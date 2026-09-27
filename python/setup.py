"""Generate the gRPC stubs from ../protos/Kraken2.proto at build time.

The proto file stays the single source of truth; nothing generated is
committed. `pip install .` or `pip install -e .` both run this step.
"""

import os
import subprocess
import sys

from setuptools import setup
from setuptools.command.build_py import build_py

HERE = os.path.abspath(os.path.dirname(__file__))
PROTO_DIR = os.path.normpath(os.path.join(HERE, "..", "protos"))
PACKAGE_DIR = os.path.join(HERE, "kraken2_client")


def generate_stubs():
    proto = os.path.join(PROTO_DIR, "Kraken2.proto")
    if not os.path.exists(proto):
        raise SystemExit("Cannot find %s; build from the repository checkout." % proto)
    cmd = [
        sys.executable, "-m", "grpc_tools.protoc",
        "-I", PROTO_DIR,
        "--python_out", PACKAGE_DIR,
        "--grpc_python_out", PACKAGE_DIR,
        proto,
    ]
    subprocess.check_call(cmd)
    # grpc_tools writes an absolute import; make it relative to the package.
    grpc_file = os.path.join(PACKAGE_DIR, "Kraken2_pb2_grpc.py")
    with open(grpc_file) as f:
        text = f.read()
    text = text.replace("import Kraken2_pb2 as", "from . import Kraken2_pb2 as")
    with open(grpc_file, "w") as f:
        f.write(text)


class BuildWithStubs(build_py):
    def run(self):
        generate_stubs()
        super().run()


setup(cmdclass={"build_py": BuildWithStubs})
