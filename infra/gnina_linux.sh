#!/bin/bash
# Run the GNINA x86-64 Linux binary natively (hosted / Linux). Usage: [GNINA_TIMEOUT=secs] infra/gnina_linux.sh <gnina args>
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
BIN="${GNINA_BIN:-$ROOT/bin/gnina.cuda12.8.static}"
exec env LD_LIBRARY_PATH="${CUDA_LIBS:-$ROOT/infra/cudalibs}:$LD_LIBRARY_PATH" ${GNINA_TIMEOUT:+timeout -k 10 $GNINA_TIMEOUT} "$BIN" "$@"
