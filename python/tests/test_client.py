"""End-to-end tests for the Python client against a real kraken2_server.

Environment variables (defaults are relative to the repository root):
  KRAKEN2_SERVER  path to kraken2_server (default build/server/kraken2_server)
  KRAKEN2_CLIENT  path to the C++ kraken2_client, used as the reference
  KRAKEN2_DB      kraken2 database (default testing/virus-zymo-kraken2)
  KRAKEN2_READS   FASTQ used for the tests (default testing/virus-zymo-kraken.reads.fastq.gz)
"""

import gzip
import io
import os
import shutil
import socket
import subprocess
import time

import pytest

from kraken2_client import (TLS, Classification, Client, ShutdownRefused, read_fastx, read_pairs)

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
SERVER = os.environ.get("KRAKEN2_SERVER", os.path.join(ROOT, "build", "server", "kraken2_server"))
CPP_CLIENT = os.environ.get("KRAKEN2_CLIENT", os.path.join(ROOT, "build", "client", "kraken2_client"))
DB = os.environ.get("KRAKEN2_DB", os.path.join(ROOT, "testing", "virus-zymo-kraken2"))
READS = os.environ.get("KRAKEN2_READS", os.path.join(ROOT, "testing", "virus-zymo-kraken.reads.fastq.gz"))

pytestmark = pytest.mark.skipif(
    not (os.path.exists(SERVER) and os.path.isdir(DB) and os.path.exists(READS)),
    reason="kraken2_server binary, test database or reads not found")


def free_port():
    s = socket.socket()
    s.bind(("127.0.0.1", 0))
    port = s.getsockname()[1]
    s.close()
    return port


class Server:
    def __init__(self, *extra_args):
        self.port = free_port()
        self.proc = subprocess.Popen(
            [SERVER, "--db", DB, "--host-ip", "127.0.0.1", "--port", str(self.port), "--thread-pool", "2"] + list(extra_args),
            stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    def stop(self):
        if self.proc.poll() is None:
            self.proc.terminate()
        return self.proc.wait(timeout=30)


@pytest.fixture(scope="module")
def server():
    s = Server("--allow-remote-shutdown")
    with Client("127.0.0.1", s.port) as c:
        c.wait_ready(poll_seconds=0.5, timeout=60)
    yield s
    s.stop()


@pytest.fixture(scope="module")
def pairs(tmp_path_factory):
    """Split the test reads into two mate files with matching names."""
    d = tmp_path_factory.mktemp("pairs")
    r1, r2 = d / "r1.fq", d / "r2.fq"
    recs = list(read_fastx(READS))
    n = min(1500, len(recs) // 2)
    with open(r1, "w") as f1, open(r2, "w") as f2:
        for i in range(n):
            a, b = recs[i], recs[n + i]
            f1.write("@frag%d/1\n%s\n+\n%s\n" % (i + 1, a[1], a[2]))
            f2.write("@frag%d/2\n%s\n+\n%s\n" % (i + 1, b[1], b[2]))
    return str(r1), str(r2)


def cpp_client(port, *args):
    return subprocess.run([CPP_CLIENT, "--host-ip", "127.0.0.1", "--port", str(port)] + list(args),
                          check=True, capture_output=True, text=True).stdout


def test_read_fastx_parses_fasta_and_fastq(tmp_path):
    fa = tmp_path / "x.fa"
    fa.write_text(">a desc\nACGT\nGG\n>b\nTT\n")
    assert list(read_fastx(str(fa))) == [("a", "ACGTGG", ""), ("b", "TT", "")]
    fq = tmp_path / "x.fq.gz"
    with gzip.open(fq, "wt") as f:
        f.write("@r1 comment\nACGT\n+\nIIII\n@r2\nGG\n+\n##\n")
    assert list(read_fastx(str(fq))) == [("r1", "ACGT", "IIII"), ("r2", "GG", "##")]


def test_read_pairs_rejects_mismatched_lengths(tmp_path):
    a, b = tmp_path / "a.fq", tmp_path / "b.fq"
    a.write_text("@x/1\nAC\n+\nII\n@y/1\nAC\n+\nII\n")
    b.write_text("@x/2\nGT\n+\nII\n")
    it = read_pairs(str(a), str(b))
    assert next(it) == ("x/1", "AC", "II", "GT", "II")
    with pytest.raises(ValueError):
        next(it)


def test_classification_line_format():
    c = Classification(True, "r", 562, 100, "562:5 0:2")
    assert c.to_kraken_line() == "C\tr\t562\t100\t562:5 0:2"
    p = Classification(False, "f", 0, 100, "0:0", paired=True, mate_length=90)
    assert p.to_kraken_line() == "U\tf\t0\t100|90\t0:0"


def test_ready_and_summary(server):
    with Client("127.0.0.1", server.port) as c:
        assert c.ready()
        assert isinstance(c.summary(), str)


def test_single_end_matches_cpp_client(server):
    with Client("127.0.0.1", server.port) as c:
        report = io.StringIO()
        lines = [h.to_kraken_line() for h in c.classify(read_fastx(READS), report=report)]
    ref = cpp_client(server.port, "--sequence", READS).rstrip("\n").split("\n")
    assert len(lines) == len(ref) > 0
    assert lines == ref
    assert report.getvalue().startswith("% of Seqs")


def test_paired_end_matches_cpp_client(server, pairs):
    r1, r2 = pairs
    with Client("127.0.0.1", server.port) as c:
        lines = [h.to_kraken_line() for h in c.classify(read_pairs(r1, r2))]
    ref = cpp_client(server.port, "--sequence", r1, "--sequence2", r2).rstrip("\n").split("\n")
    assert lines == ref
    assert all("|:|" in line for line in lines)
    assert lines[0].split("\t")[1] == "frag1"


def test_results_arrive_in_input_order(server):
    ids = [r[0] for r in read_fastx(READS)]
    with Client("127.0.0.1", server.port) as c:
        got = [h.read_id for h in c.classify(read_fastx(READS))]
    assert got == ids


def test_empty_input(server):
    with Client("127.0.0.1", server.port) as c:
        assert list(c.classify([])) == []


def test_shutdown_refused_without_flag():
    s = Server()
    try:
        with Client("127.0.0.1", s.port) as c:
            c.wait_ready(poll_seconds=0.5, timeout=60)
            with pytest.raises(ShutdownRefused):
                c.shutdown()
        assert s.proc.poll() is None
    finally:
        assert s.stop() == 0


@pytest.mark.skipif(shutil.which("openssl") is None, reason="openssl not available")
def test_tls(tmp_path):
    crt, key = tmp_path / "s.crt", tmp_path / "s.key"
    subprocess.run(["openssl", "req", "-x509", "-newkey", "rsa:2048", "-nodes", "-days", "2",
                    "-subj", "/CN=localhost", "-addext", "subjectAltName=DNS:localhost,IP:127.0.0.1",
                    "-keyout", str(key), "-out", str(crt)], check=True, capture_output=True)
    s = Server("--tls-cert", str(crt), "--tls-key", str(key), "--allow-remote-shutdown")
    try:
        with Client("127.0.0.1", s.port, TLS(ca=str(crt))) as c:
            c.wait_ready(poll_seconds=0.5, timeout=60)
            n = sum(1 for _ in c.classify(read_fastx(READS)))
            assert n > 0
            assert c.shutdown()
        assert s.proc.wait(timeout=30) == 0
    finally:
        s.stop()


def test_remote_shutdown(server):
    # last test in the module: stops the shared server
    with Client("127.0.0.1", server.port) as c:
        assert c.shutdown()
    assert server.proc.wait(timeout=30) == 0
