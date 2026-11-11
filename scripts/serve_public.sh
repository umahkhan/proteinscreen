#!/bin/bash
# Serve pscreen to the public from THIS Mac through a free Cloudflare quick tunnel. Ctrl-C stops everything.
#   scripts/serve_public.sh              start server + tunnel, print the public link
#   scripts/serve_public.sh --no-tunnel  start the hosted-mode server on this Mac only (127.0.0.1) for testing
# Needs once:  brew install cloudflared
# Visitors get private workspaces, size/run caps and a queue (hosted mode). Their uploads are processed on THIS computer: see docs/SERVE_FROM_MAC.md.
set -u
cd "$(dirname "$0")/.."
PORT="${PORT:-8601}"; MODE="${1:-tunnel}"
export PSCREEN_HOSTED=1 PSCREEN_DATA="$PWD/.public_data" PSCREEN_LIMITS="${PSCREEN_LIMITS:-off}" PSCREEN_SLOTS="${PSCREEN_SLOTS:-1}" PSCREEN_RETENTION_H="${PSCREEN_RETENTION_H:-48}" \
       PSCREEN_MAX_COMPOUNDS="${PSCREEN_MAX_COMPOUNDS:-25}" PSCREEN_MAX_QUEUE="${PSCREEN_MAX_QUEUE:-6}"
mkdir -p "$PSCREEN_DATA"
PIDS=()
cleanup() { echo; echo "Stopping..."; for p in "${PIDS[@]}"; do kill "$p" 2>/dev/null; done; pkill -f "streamlit run app.py --server.port $PORT" 2>/dev/null; exit 0; }
trap cleanup INT TERM

limactl start gnina >/dev/null 2>&1 || true                       # GNINA runs inside this VM on macOS
caffeinate -dimsu -w $$ &                                         # keep the Mac awake while this script runs
PIDS+=($!)

if [ "$MODE" != "--no-tunnel" ]; then
  command -v cloudflared >/dev/null || { echo "Install the tunnel tool once:  brew install cloudflared"; exit 1; }
  : > "$PSCREEN_DATA/cloudflared.log"
  cloudflared tunnel --url "http://127.0.0.1:$PORT" --no-autoupdate > "$PSCREEN_DATA/cloudflared.log" 2>&1 &
  PIDS+=($!)
  URL=""
  for _ in $(seq 1 40); do URL=$(grep -o 'https://[a-z0-9-]*\.trycloudflare\.com' "$PSCREEN_DATA/cloudflared.log" | head -1); [ -n "$URL" ] && break; sleep 1; done
  [ -n "$URL" ] || { echo "Could not get a public link. See $PSCREEN_DATA/cloudflared.log"; cleanup; }
  export PSCREEN_PUBLIC_URL="$URL"
fi

pkill -f "streamlit run app.py --server.port $PORT" 2>/dev/null; sleep 1
streamlit run app.py --server.port "$PORT" --server.address 127.0.0.1 --server.headless true --server.enableXsrfProtection true \
  --server.enableCORS true --server.maxUploadSize 500 --browser.gatherUsageStats false > "$PSCREEN_DATA/streamlit.log" 2>&1 &
PIDS+=($!)
for _ in $(seq 1 30); do curl -sf "http://127.0.0.1:$PORT/_stcore/health" >/dev/null && break; sleep 1; done
echo
echo "pscreen is running."
echo "  This Mac only : http://127.0.0.1:$PORT"
[ "$MODE" != "--no-tunnel" ] && echo "  PUBLIC LINK   : $URL     (share this; it changes every time you restart)"
echo "Keep this window open and the Mac awake. Press Ctrl-C to stop."
wait
