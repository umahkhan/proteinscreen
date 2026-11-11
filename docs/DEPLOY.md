# Hosting pscreen for other people

Hosted mode (`PSCREEN_HOSTED=1`) gives every visitor a private workspace (`?ws=<token>` in the page address), a read-only example run, and automatic
deletion of idle workspaces (`PSCREEN_RETENTION_H`, default 48 h). Local use is unchanged.

## Option 1: serve from your own computer (used so far)
`scripts/serve_public.sh` starts the full node-based app in hosted mode behind a free Cloudflare quick tunnel. See `docs/SERVE_FROM_MAC.md`.
By default it sets `PSCREEN_LIMITS=off` (no caps, no queue). Set `PSCREEN_LIMITS=on` to apply the free-tier caps in `pscreen/hosting.py`
(25 compounds, 2 covalent residues, 1 active run per workspace, a global queue and a 10-minute job timeout).

## Option 2: a container on a cloud host (not verified)
The `Dockerfile` builds the app with GNINA and the tool environments. It has never been built or run. Hugging Face Docker and Gradio Spaces need a
paid plan; static Spaces are free but cannot run Python. Pocket finding alone can run on a free Streamlit host: see `pockets_app.py` and `deploy/cloud/make_repo.sh`.

## Not yet verified
- The Docker image has not been built (no Docker on the development Mac).
- Shared hosts with 2 CPUs are several times slower than the Mac used so far.
- Licences: GNINA Apache-2.0, P2Rank MIT, Lacuna MIT, ADMET-AI MIT, PocketMiner MIT per its README (no LICENSE file detected by GitHub: confirm before launch).
- Uploads can contain proprietary structures. The app states that the operator can read stored files.
