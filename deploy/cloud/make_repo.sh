#!/bin/bash
# Assemble the minimal repository for Streamlit Community Cloud (pockets-only app) into ~/pscreen-pockets-repo.
set -e
SRC="$(cd "$(dirname "$0")/../.." && pwd)"; OUT="${1:-$HOME/pscreen-pockets-repo}"
rm -rf "$OUT"; mkdir -p "$OUT/.streamlit"
cp "$SRC/pockets_app.py" "$OUT/streamlit_app.py"
cp -r "$SRC/pscreen" "$SRC/static" "$SRC/examples_data" "$OUT/"
rm -rf "$OUT/pscreen/__pycache__" "$OUT"/pscreen/*/__pycache__
# PocketMiner (MIT): source + the original checkpoint only (19 MB; the 500 MB training data is not needed)
PM="$SRC/third_party/ae-pocketminer"; mkdir -p "$OUT/third_party/ae-pocketminer/models"
cp -r "$PM/src" "$OUT/third_party/ae-pocketminer/src"; rm -rf "$OUT/third_party/ae-pocketminer/src/__pycache__" "$OUT/third_party/ae-pocketminer/src/training-yaml-files"
cp "$PM"/models/pocketminer.* "$OUT/third_party/ae-pocketminer/models/"; cp "$PM/LICENSE" "$OUT/third_party/ae-pocketminer/"
cp "$SRC/deploy/cloud/requirements.txt" "$SRC/deploy/cloud/packages.txt" "$OUT/"
# Cloud-safe config: theme and bundled font only (no local-only server settings)
python3 - "$SRC/.streamlit/config.toml" "$OUT/.streamlit/config.toml" <<'PY'
import sys, re
t = open(sys.argv[1]).read(); t = re.sub(r"\[server\].*?(?=\n\[|\Z)", "[server]\nenableStaticServing = true\n", t, flags=re.S); open(sys.argv[2], "w").write(t)
PY
cat > "$OUT/README.md" <<'MD'
# pscreen · Pockets
Find candidate binding pockets (conventional and cryptic) on a protein structure. Pocket finding only, no docking.
Uses P2Rank, Lacuna and (when installed) PocketMiner. Predictions only: not evidence of binding or druggability.
MD
echo "assembled $OUT"; du -sh "$OUT"
