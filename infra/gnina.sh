#!/bin/bash
# Run GNINA (x86-64 Linux binary) inside the Lima VM via Rosetta.
# Usage: [GNINA_TIMEOUT=secs] infra/gnina.sh <gnina args>   (all paths must be under $HOME/proteinscreen)
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
limactl start gnina >/dev/null 2>&1 || true
exec limactl shell --workdir "$ROOT" gnina -- env LD_LIBRARY_PATH="$ROOT/infra/cudalibs" \
  ${GNINA_TIMEOUT:+timeout -k 10 $GNINA_TIMEOUT} "$ROOT/bin/gnina.cuda12.8.static" "$@"
