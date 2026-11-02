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

import collections
import gzip
import io
import os
import sys
import threading
import time
import zlib
from dataclasses import dataclass
from typing import Callable, Iterable, Iterator, Optional, Tuple, Union

import grpc

from . import Kraken2_pb2 as pb
from . import Kraken2_pb2_grpc as pb_grpc

__all__ = ["Client", "Classification", "Record", "PairedRecord", "read_fastx", "read_pairs",
           "TLS", "InputError", "ServerNotReady", "ShutdownRefused"]

Record = Tuple[str, str, str]
PairedRecord = Tuple[str, str, str, str, str]

# Buffering is bounded in bytes as well as records so that long reads do not
# make the client hold its whole input in memory.
BATCH_SIZE = 2000                        # at most this many records per gRPC message
BATCH_BYTES = 8 * 1024 * 1024            # or this much sequence and quality data
MAX_IN_FLIGHT = 64000                    # records sent but not yet answered
MAX_IN_FLIGHT_BYTES = 256 * 1024 * 1024  # bytes sent but not yet answered
MAX_MESSAGE_BYTES = 128 * 1024 * 1024    # the server's message limit
FLUSH_SECONDS = 1.0     # send a partial batch once its first record is this old
FEED_CHUNK = 512                        # records handed between threads at once
FEED_CHUNK_BYTES = 8 * 1024 * 1024      # or this much sequence
FEED_AHEAD = 8                          # chunks read ahead of the request generator
KEEPALIVE_TIMEOUT_MS = 20000   # a keepalive ping unanswered this long ends the stream

_END = object()


class ServerNotReady(RuntimeError):
    """The server is up but its database is not loaded yet."""


class ShutdownRefused(PermissionError):
    """The server was not started with --allow-remote-shutdown."""


class InputError(ValueError):
    """An input file could not be parsed, for example a FASTQ file that
    ends part way through a record or a truncated gzip file."""


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


class _InFlight:
    """Records and bytes sent but not yet answered. The request generator
    acquires room for each record (blocking on gRPC's thread while the
    limits are exceeded) and the result loop releases records in order.
    close() unblocks the generator when classify() exits early."""

    def __init__(self):
        self._cond = threading.Condition()
        self._sizes = collections.deque()   # bytes of each record in flight, oldest first
        self._bytes = 0
        self._closed = False

    def acquire(self, nbytes: int):
        with self._cond:
            while not self._closed and self._sizes and (
                    len(self._sizes) >= MAX_IN_FLIGHT or self._bytes + nbytes > MAX_IN_FLIGHT_BYTES):
                self._cond.wait(1.0)
            self._sizes.append(nbytes)
            self._bytes += nbytes

    def release(self):
        with self._cond:
            if self._sizes:
                self._bytes -= self._sizes.popleft()
            self._cond.notify_all()

    def close(self):
        with self._cond:
            self._closed = True
            self._cond.notify_all()


class _Feed:
    """Records handed from the thread that reads the caller's iterable to
    gRPC's request generator, in chunks so that the two threads meet once
    per FEED_CHUNK records rather than once per record. The reader adds
    records to a pending chunk and queues it when full; the generator takes
    queued chunks, and the pending chunk too once nothing has been queued
    for its timeout, so records that arrive slowly are not held back.
    At most FEED_AHEAD chunks wait to be taken."""

    def __init__(self):
        self.lock = threading.Lock()
        self.pending = []
        self.pending_bytes = 0
        self.chunks = collections.deque()
        self.ready = threading.Event()   # a chunk was queued for a waiting generator
        self.space = threading.Event()   # a chunk was taken for a waiting reader
        self.consumer_waiting = False
        self.producer_waiting = False
        self.stop = threading.Event()

    def put(self, item, nbytes: int) -> bool:
        """Add one item. False once the feed has been stopped."""
        with self.lock:
            self.pending.append(item)
            self.pending_bytes += nbytes
            if len(self.pending) < FEED_CHUNK and self.pending_bytes < FEED_CHUNK_BYTES:
                return not self.stop.is_set()
            chunk, self.pending, self.pending_bytes = self.pending, [], 0
        return self._queue(chunk)

    def close(self):
        """Queue whatever is pending; the reader calls this at the end."""
        with self.lock:
            chunk, self.pending, self.pending_bytes = self.pending, [], 0
        if chunk:
            self._queue(chunk)

    def _queue(self, chunk) -> bool:
        while len(self.chunks) >= FEED_AHEAD and not self.stop.is_set():
            self.producer_waiting = True
            self.space.clear()
            if len(self.chunks) >= FEED_AHEAD:
                self.space.wait(0.5)
            self.producer_waiting = False
        if self.stop.is_set():
            return False
        self.chunks.append(chunk)
        if self.consumer_waiting:
            self.ready.set()
        return True

    def get(self, timeout: float):
        """The next chunk of items, or None if nothing arrived in time."""
        if not self.chunks:
            self.consumer_waiting = True
            self.ready.clear()
            if not self.chunks:
                self.ready.wait(timeout)
            self.consumer_waiting = False
        with self.lock:
            # A queued chunk is always older than the pending one: the
            # reader queues a full chunk before it adds the next item.
            if self.chunks:
                chunk = self.chunks.popleft()
            elif self.pending:
                chunk, self.pending, self.pending_bytes = self.pending, [], 0
            else:
                return None
        if self.producer_waiting:
            self.space.set()
        return chunk


class Client:
    def __init__(self, host: str = "localhost", port: int = 8080, tls: Optional[TLS] = None,
                 keepalive: int = 60):
        """keepalive: seconds without data on a stream before the client
        pings the server; a ping unanswered for 20 s ends the stream with
        UNAVAILABLE, so a frozen server or a network that drops without a
        reset is noticed. At least 10 (the server refuses pings more often);
        0 disables."""
        if keepalive and keepalive < 10:
            raise ValueError("keepalive must be 0 or at least 10 seconds")
        target = "%s:%d" % (host, port)
        options = [("grpc.max_receive_message_length", -1),
                   ("grpc.max_send_message_length", MAX_MESSAGE_BYTES)]
        if keepalive:
            options += [("grpc.keepalive_time_ms", keepalive * 1000),
                        ("grpc.keepalive_timeout_ms", KEEPALIVE_TIMEOUT_MS),
                        # current gRPC times keepalive pings out with this
                        ("grpc.http2.ping_timeout_ms", KEEPALIVE_TIMEOUT_MS),
                        ("grpc.keepalive_permit_without_calls", 0),
                        ("grpc.http2.max_pings_without_data", 0)]
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
                 report: Optional[io.TextIOBase] = None,
                 on_message: Optional[Callable[[], None]] = None) -> Iterator[Classification]:
        """Stream records to the server and yield one Classification per
        record, in input order. If `report` is a writable text stream, the
        kraken2 style report for this stream is written to it at the end.

        Records are sent in messages of at most BATCH_SIZE records or
        BATCH_BYTES of sequence, so long reads never exceed the server's
        message limit, and at most MAX_IN_FLIGHT records or
        MAX_IN_FLIGHT_BYTES are awaiting results, which bounds memory use
        whatever the read length. A single record over the message limit
        raises ValueError. An iterable that may block (stdin, a pipe, any
        iterable other than read_fastx or read_pairs on regular files) is
        read on a separate thread, and a partial message is sent once its
        first record is FLUSH_SECONDS old, so records that arrive slowly
        (from a running sequencer) are classified as they come rather than
        once a message is full.
        on_message, if given, is called after the results of each message
        from the server have been yielded; the command line client uses it
        to flush its output, so results reach a pipe as they arrive.
        An exception raised by the `records` iterable (for example from
        read_pairs on mismatched inputs) is re-raised after the results for
        the records already sent have been yielded."""
        in_flight = _InFlight()
        summary = []
        # An exception from the records iterable. It cannot propagate out
        # of the request generator, which gRPC runs on its own thread and
        # would report as an opaque UNKNOWN status, so it is kept here and
        # raised once the stream has been drained.
        failure = []

        # An iterable that may block (stdin, a pipe, a generator fed by a
        # running sequencer) is read on its own thread, so that the request
        # generator can send a partial message while it blocks. Records
        # read from a regular file by read_fastx or read_pairs cannot
        # block, and are read directly, which is faster.
        feed = _Feed()
        stop = feed.stop
        may_block = getattr(records, "may_block", True)

        def feeder():
            try:
                for rec in records:
                    if not feed.put(rec, _record_bytes(rec)):
                        return
            except Exception as e:  # noqa: BLE001 - re-raised by classify()
                failure.append(e)
            feed.put(_END, 0)
            feed.close()

        if may_block:
            threading.Thread(target=feeder, name="kraken2-client-feeder", daemon=True).start()

        def records_from_feed():
            # yields records, or None when nothing arrived before the
            # current batch is due (the generator then sends it)
            while True:
                if len(state["batch"].seqs):
                    timeout = max(0.0, state["started"] + FLUSH_SECONDS - time.monotonic())
                else:
                    timeout = 0.5
                chunk = feed.get(timeout)
                if chunk is None:
                    if stop.is_set():
                        return
                    yield None
                    continue
                for rec in chunk:
                    if rec is _END:
                        return
                    yield rec

        def records_direct():
            try:
                yield from records
            except Exception as e:  # noqa: BLE001 - re-raised by classify()
                failure.append(e)

        state = {"batch": pb.Kraken2SequenceRequestMulti(), "started": 0.0}

        def requests():
            batch = state["batch"]
            batch_bytes = 0
            for rec in (records_from_feed() if may_block else records_direct()):
                if rec is None:
                    if len(batch.seqs):
                        yield batch
                        batch = state["batch"] = pb.Kraken2SequenceRequestMulti()
                        batch_bytes = 0
                    continue
                rec_bytes = _record_bytes(rec)
                if rec_bytes > MAX_MESSAGE_BYTES:
                    failure.append(ValueError("record %s is larger than the %d MB message limit"
                                              % (rec[0], MAX_MESSAGE_BYTES // (1024 * 1024))))
                    break
                if len(batch.seqs) and (len(batch.seqs) >= BATCH_SIZE
                                        or batch_bytes + rec_bytes > BATCH_BYTES):
                    yield batch
                    batch = state["batch"] = pb.Kraken2SequenceRequestMulti()
                    batch_bytes = 0
                in_flight.acquire(rec_bytes)
                if not len(batch.seqs):
                    state["started"] = time.monotonic()
                r = batch.seqs.add()
                _fill(r, rec[0], rec[1], rec[2])
                if len(rec) == 5:
                    _fill(r.mate, rec[0], rec[3], rec[4])
                batch_bytes += rec_bytes
            if len(batch.seqs):
                yield batch

        try:
            for result in self._stub.ClassifyStream(requests()):
                if result.HasField("classifications"):
                    for c in result.classifications.classes:
                        in_flight.release()
                        yield Classification(c.classified, c.id.decode("utf-8", "surrogateescape"),
                                             c.tax_id, c.size, c.hitlist, c.paired, c.size2)
                    if on_message is not None:
                        on_message()
                elif result.HasField("summary"):
                    summary.append(result.summary)
        finally:
            stop.set()
            in_flight.close()
        # The report covers the records that were classified, also when
        # the record source failed part way.
        if report is not None:
            report.write("".join(summary))
        if failure:
            raise failure[0]


def _record_bytes(rec) -> int:
    """Approximate serialized size of a record: the strings plus a small
    allowance for field tags and lengths."""
    n = len(rec[0]) + len(rec[1]) + len(rec[2]) + 16
    if len(rec) == 5:
        n += len(rec[0]) + len(rec[3]) + len(rec[4]) + 16
    return n


# Read ids, sequences and qualities travel as bytes. Text read from files
# is decoded with surrogateescape, so a byte that is not UTF-8 is carried
# through unchanged and comes back the same in the result's read_id.
def _fill(msg, read_id: str, seq: str, quals: str):
    msg.id = read_id.encode("utf-8", "surrogateescape")
    msg.seq = seq.encode("utf-8", "surrogateescape")
    if quals:
        msg.format = pb.Kraken2SequenceRequest.FORMAT_FASTQ
        msg.quals = quals.encode("utf-8", "surrogateescape")
    else:
        msg.format = pb.Kraken2SequenceRequest.FORMAT_FASTA


def _open_text(path: str):
    """A text stream over a FASTA/FASTQ file, plain or gzipped, and a
    function that releases what was opened here. Works on inputs that
    cannot seek (stdin, a named pipe, process substitution): the gzip magic
    is looked at with peek(), which does not consume it. Standard input
    is left open."""
    if path == "-" and not hasattr(sys.stdin, "buffer"):
        # stdin replaced by a text stream (a test, a notebook)
        return sys.stdin, lambda: None
    raw = sys.stdin.buffer if path == "-" else open(path, "rb")
    head = raw.peek(2)[:2] if hasattr(raw, "peek") else b""
    binary = gzip.GzipFile(fileobj=raw, mode="rb") if head == b"\x1f\x8b" else raw
    text = io.TextIOWrapper(binary, encoding="utf-8", errors="surrogateescape")

    def release():
        if path == "-":
            text.detach()   # closing the wrapper would close stdin
            if binary is not raw:
                binary.close()   # a GzipFile does not close the stream it wraps
        else:
            text.close()
            raw.close()
    return text, release


class _Records:
    """An iterator of records that also says whether reading it may block
    for a long time (stdin or a pipe), as opposed to a regular file."""

    def __init__(self, it, may_block: bool):
        self._it = it
        self.may_block = may_block

    def __iter__(self):
        return self

    def __next__(self):
        return next(self._it)


def read_fastx(path: str) -> Iterator[Record]:
    """Return an iterator of (id, seq, quals) from a FASTA or FASTQ file,
    plain or gzipped, or from stdin when path is "-". The id is the header
    up to the first whitespace. quals is empty for FASTA.

    The file is opened at once, so a missing file raises OSError here
    rather than when the records are first read. A FASTQ record cut short
    (a truncated file), a truncated or corrupt gzip file, or a read error
    raises InputError when it is reached. Named pipes and gzip on stdin
    work. Bytes that are not UTF-8 are kept (decoded with surrogateescape).
    Standard input is left open."""
    f, release = _open_text(path)
    may_block = path == "-" or not os.path.isfile(path)
    return _Records(_fastx_records(f, path, release), may_block)


def _fastx_records(f, path: str, release) -> Iterator[Record]:
    try:
        header = None
        seq_lines = []
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            if line[0] == "@" and header is None:
                read_id = line[1:].split(None, 1)[0]
                seq = f.readline().rstrip("\n")
                plus = f.readline()
                quals = f.readline().rstrip("\n")
                if not plus.startswith("+") or len(quals) != len(seq):
                    raise InputError(
                        "%s: record %s is incomplete; the file appears to be truncated" % (path, read_id))
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
    except (EOFError, OSError, zlib.error) as e:
        # gzip raises EOFError for a compressed stream that ends early,
        # BadGzipFile (an OSError) for a bad CRC or header, and zlib.error
        # for corrupt data; a read error is an OSError too.
        raise InputError("%s: %s" % (path, e)) from None
    finally:
        release()


def read_pairs(path1: str, path2: str) -> Iterator[PairedRecord]:
    """Return an iterator of paired records from two files matched by
    position. Both files are opened at once. Stops at the end of the
    shorter file and raises ValueError if the lengths differ."""
    it1 = read_fastx(path1)
    it2 = read_fastx(path2)
    return _Records(_pairs(it1, it2, path1, path2), it1.may_block or it2.may_block)


def _pairs(it1, it2, path1: str, path2: str) -> Iterator[PairedRecord]:
    while True:
        r1 = next(it1, None)
        r2 = next(it2, None)
        if r1 is None and r2 is None:
            return
        if r1 is None or r2 is None:
            raise ValueError("paired-end inputs have different numbers of reads: %s, %s" % (path1, path2))
        yield (r1[0], r1[1], r1[2], r2[1], r2[2])
