"""Python client for kraken2_server.

Example:

    from kraken2_client import Client, read_fastx

    with Client("localhost", 8080) as client:
        for hit in client.classify(read_fastx("reads.fq.gz")):
            print(hit.classified, hit.read_id, hit.tax_id, hit.hitlist)

Records are (id, seq, quals) tuples, or (id, seq, quals, mate_seq, mate_quals)
for read pairs; quals is an empty string for FASTA. Results are yielded in the
order the records were sent.
"""

import gzip
import io
import threading
from dataclasses import dataclass
from typing import Iterable, Iterator, Optional, Tuple, Union

import grpc

from . import Kraken2_pb2 as pb
from . import Kraken2_pb2_grpc as pb_grpc

__all__ = ["Client", "Classification", "Record", "PairedRecord", "read_fastx", "TLS",
           "ServerNotReady", "ShutdownRefused"]

Record = Tuple[str, str, str]
PairedRecord = Tuple[str, str, str, str, str]

BATCH_SIZE = 2000        # records per gRPC message, as in the C++ client
MAX_IN_FLIGHT = 64000    # records sent but not yet answered
MAX_MESSAGE_BYTES = 128 * 1024 * 1024


class ServerNotReady(RuntimeError):
    """The server is up but its database is not loaded yet."""


class ShutdownRefused(PermissionError):
    """The server was not started with --allow-remote-shutdown."""


@dataclass
class TLS:
    """TLS settings. `ca` is a PEM CA bundle (or the server's self-signed
    certificate); None uses the system roots. `cert` and `key` present a
    client certificate for mutual TLS. `server_name` overrides the name
    checked in the server certificate."""
    ca: Optional[str] = None
    cert: Optional[str] = None
    key: Optional[str] = None
    server_name: Optional[str] = None


@dataclass
class Classification:
    classified: bool
    read_id: str
    tax_id: int
    length: int
    hitlist: str
    paired: bool = False
    mate_length: int = 0

    def to_kraken_line(self) -> str:
        """The record formatted as one line of kraken2 output."""
        length = "%d|%d" % (self.length, self.mate_length) if self.paired else str(self.length)
        return "\t".join(["C" if self.classified else "U", self.read_id,
                          str(self.tax_id), length, self.hitlist])


def _read_pem(path: Optional[str]) -> Optional[bytes]:
    if path is None:
        return None
    with open(path, "rb") as f:
        return f.read()


class Client:
    def __init__(self, host: str = "localhost", port: int = 8080, tls: Optional[TLS] = None):
        target = "%s:%d" % (host, port)
        options = [("grpc.max_receive_message_length", -1),
                   ("grpc.max_send_message_length", MAX_MESSAGE_BYTES)]
        if tls is None:
            self._channel = grpc.insecure_channel(target, options=options)
        else:
            creds = grpc.ssl_channel_credentials(
                root_certificates=_read_pem(tls.ca),
                private_key=_read_pem(tls.key),
                certificate_chain=_read_pem(tls.cert))
            if tls.server_name:
                options.append(("grpc.ssl_target_name_override", tls.server_name))
            self._channel = grpc.secure_channel(target, creds, options=options)
        self._stub = pb_grpc.Kraken2ServiceStub(self._channel)

    def close(self):
        self._channel.close()

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()

    def ready(self) -> bool:
        """True when the server has loaded its database. Raises grpc.RpcError
        for transport failures."""
        try:
            return self._stub.ServerReady(pb.Kraken2ReadyRequest()).ready
        except grpc.RpcError as e:
            if e.code() == grpc.StatusCode.UNAVAILABLE and e.details().startswith("Index not loaded"):
                return False
            raise

    def wait_ready(self, poll_seconds: float = 2.0, timeout: Optional[float] = None,
                   connect_timeout: float = 30.0):
        """Block until ready(). Transport failures (connection refused, TLS
        mismatch) are retried for connect_timeout seconds and then raised as
        grpc.RpcError. A server whose database is still loading is waited for
        until timeout seconds, then ServerNotReady is raised; None waits
        indefinitely."""
        import time
        start = time.monotonic()
        deadline = None if timeout is None else start + timeout
        while True:
            try:
                if self.ready():
                    return
            except grpc.RpcError as e:
                if e.code() != grpc.StatusCode.UNAVAILABLE or time.monotonic() - start > connect_timeout:
                    raise
            if deadline is not None and time.monotonic() > deadline:
                raise ServerNotReady("server database not loaded after %.0f s" % timeout)
            time.sleep(poll_seconds)

    def summary(self) -> str:
        """The server's cumulative kraken2 style report, or a notice if the
        server was started with --no-stats."""
        return self._stub.GetSummary(pb.Kraken2SummaryRequest()).summary

    def shutdown(self) -> bool:
        try:
            return self._stub.RemoteShutdown(pb.Kraken2ShutdownRequest()).successful
        except grpc.RpcError as e:
            if e.code() == grpc.StatusCode.PERMISSION_DENIED:
                raise ShutdownRefused(e.details()) from None
            raise

    def classify(self, records: Iterable[Union[Record, PairedRecord]],
                 report: Optional[io.TextIOBase] = None) -> Iterator[Classification]:
        """Stream records to the server and yield one Classification per
        record, in input order. If `report` is a writable text stream, the
        kraken2 style report for this stream is written to it at the end."""
        in_flight = threading.Semaphore(MAX_IN_FLIGHT)
        summary = []

        def requests():
            batch = pb.Kraken2SequenceRequestMulti()
            for rec in records:
                in_flight.acquire()
                r = batch.seqs.add()
                _fill(r, rec[0], rec[1], rec[2])
                if len(rec) == 5:
                    _fill(r.mate, rec[0], rec[3], rec[4])
                if len(batch.seqs) >= BATCH_SIZE:
                    yield batch
                    batch = pb.Kraken2SequenceRequestMulti()
            if len(batch.seqs):
                yield batch

        for result in self._stub.ClassifyStream(requests()):
            if result.HasField("classifications"):
                for c in result.classifications.classes:
                    in_flight.release()
                    yield Classification(c.classified, c.id, c.tax_id, c.size, c.hitlist,
                                         c.paired, c.size2)
            elif result.HasField("summary"):
                summary.append(result.summary)
        if report is not None:
            report.write("".join(summary))


def _fill(msg, read_id: str, seq: str, quals: str):
    msg.id = read_id
    msg.seq = seq
    if quals:
        msg.format = pb.Kraken2SequenceRequest.FORMAT_FASTQ
        msg.quals = quals
    else:
        msg.format = pb.Kraken2SequenceRequest.FORMAT_FASTA


def _open_text(path: str):
    if path == "-":
        import sys
        return sys.stdin
    f = open(path, "rb")
    if f.read(2) == b"\x1f\x8b":
        f.seek(0)
        return io.TextIOWrapper(gzip.GzipFile(fileobj=f))
    f.seek(0)
    return io.TextIOWrapper(f)


def read_fastx(path: str) -> Iterator[Record]:
    """Yield (id, seq, quals) from a FASTA or FASTQ file, plain or gzipped,
    or from stdin when path is "-". The id is the header up to the first
    whitespace. quals is empty for FASTA."""
    with _open_text(path) as f:
        header = None
        seq_lines = []
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            if line[0] == "@" and header is None:
                read_id = line[1:].split(None, 1)[0]
                seq = next(f).rstrip("\n")
                next(f)  # '+'
                quals = next(f).rstrip("\n")
                yield (read_id, seq, quals)
            elif line[0] == ">":
                if header is not None:
                    yield (header, "".join(seq_lines), "")
                header = line[1:].split(None, 1)[0]
                seq_lines = []
            else:
                seq_lines.append(line)
        if header is not None:
            yield (header, "".join(seq_lines), "")


def read_pairs(path1: str, path2: str) -> Iterator[PairedRecord]:
    """Yield paired records from two files matched by position. Stops at the
    end of the shorter file and raises ValueError if the lengths differ."""
    it1, it2 = read_fastx(path1), read_fastx(path2)
    while True:
        r1 = next(it1, None)
        r2 = next(it2, None)
        if r1 is None and r2 is None:
            return
        if r1 is None or r2 is None:
            raise ValueError("paired-end inputs have different numbers of reads: %s, %s" % (path1, path2))
        yield (r1[0], r1[1], r1[2], r2[1], r2[2])
