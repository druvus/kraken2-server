"""Command line client with the same options as kraken2_client."""

import argparse
import sys

import grpc

from . import TLS, Client, ShutdownRefused, read_fastx, read_pairs

EX_USAGE, EX_DATAERR, EX_UNAVAILABLE = 64, 65, 69


def main(argv=None):
    ap = argparse.ArgumentParser(prog="k2client", description="Classify reads with a kraken2_server.")
    ap.add_argument("-s", "--sequence", help="FASTA/FASTQ file, plain or gzipped, or - for stdin. Omit to print the server summary.")
    ap.add_argument("-2", "--sequence2", help="mate file for paired-end reads, in the same order as --sequence")
    ap.add_argument("-r", "--report", help="write the kraken2 style report to this file")
    ap.add_argument("-i", "--host-ip", default="localhost", help="server address (default localhost)")
    ap.add_argument("-p", "--port", type=int, default=8080, help="server port (default 8080)")
    ap.add_argument("-k", "--shutdown", action="store_true", help="ask the server to shut down")
    ap.add_argument("-t", "--tls", action="store_true", help="connect with TLS using the system CA roots")
    ap.add_argument("--tls-ca", help="PEM CA bundle to verify the server (implies --tls)")
    ap.add_argument("--tls-cert", help="PEM client certificate for mutual TLS (implies --tls)")
    ap.add_argument("--tls-key", help="PEM client key for --tls-cert")
    ap.add_argument("--tls-server-name", help="name to verify in the server certificate")
    ap.add_argument("--wait", type=float, default=None, help="seconds to wait for the server database to load (default: forever)")
    args = ap.parse_args(argv)

    if args.sequence2 and not args.sequence:
        ap.error("--sequence2 requires --sequence")
    if bool(args.tls_cert) != bool(args.tls_key):
        ap.error("--tls-cert and --tls-key must be given together")

    tls = None
    if args.tls or args.tls_ca or args.tls_cert or args.tls_server_name:
        tls = TLS(ca=args.tls_ca, cert=args.tls_cert, key=args.tls_key, server_name=args.tls_server_name)

    try:
        with Client(args.host_ip, args.port, tls) as client:
            if args.shutdown:
                try:
                    ok = client.shutdown()
                except ShutdownRefused as e:
                    print("Shutdown refused: %s" % e, file=sys.stderr)
                    return 7
                print("Shutdown request %s." % ("processed" if ok else "not processed"), file=sys.stderr)
                return 0
            if not args.sequence:
                print(client.summary())
                return 0
            client.wait_ready(timeout=args.wait)
            if args.sequence2:
                records = read_pairs(args.sequence, args.sequence2)
            else:
                records = read_fastx(args.sequence)
            report = open(args.report, "w") if args.report else None
            try:
                out = sys.stdout
                for hit in client.classify(records, report=report):
                    out.write(hit.to_kraken_line())
                    out.write("\n")
            finally:
                if report:
                    report.close()
            return 0
    except ValueError as e:
        print(str(e), file=sys.stderr)
        return EX_DATAERR
    except grpc.RpcError as e:
        print("RPC failed: %s: %s" % (e.code().name, e.details()), file=sys.stderr)
        return EX_UNAVAILABLE


if __name__ == "__main__":
    sys.exit(main())
