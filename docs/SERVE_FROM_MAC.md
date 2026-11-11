# Serving pscreen from your own Mac (free)

Your Mac is the server; a free Cloudflare quick tunnel gives it a public `https://….trycloudflare.com` link. Visitors need nothing installed.

## Start / stop
1. Once: `brew install cloudflared`
2. `cd ~/proteinscreen && scripts/serve_public.sh`  → prints the **PUBLIC LINK**. Share it. Keep the window open and the Mac on power.
3. Ctrl-C stops the server, the tunnel and the keep-awake. `scripts/serve_public.sh --no-tunnel` runs the same hosted mode on this Mac only (http://127.0.0.1:8601).

## What is protected (tested)
- Each visitor gets a private random workspace (`?ws=`); forged links are replaced; one visitor cannot list another's runs.
- Upload file names are sanitised (a name like `../../x.pdb` stays inside the visitor's folder).
- No caps by default (`PSCREEN_LIMITS=off`): any number of runs may go at once and each uses the settings the visitor chooses, as in the local app. Start with `PSCREEN_LIMITS=on` to apply the free-tier caps in `pscreen/hosting.py` (25 compounds, 2 covalent residues, 1 active run per visitor, a server-wide queue, 600 s per docking job).
- Workspaces are deleted after 48 h idle (`PSCREEN_RETENTION_H`). XSRF protection is on and uploads are limited to 500 MB.

## What is NOT protected: read before sharing widely
- Uploads are parsed by RDKit, Biopython, P2Rank (Java), TensorFlow and GNINA **on your Mac under your user account**. A bug in any of them could be abused. For a wider audience use a dedicated macOS user, or another computer.
- You (the operator) can read every visitor's files. The app says so. Do not promise confidentiality.
- No per-person rate limit (a visitor can open many workspaces; with no caps the load is bounded only by your Mac). The tunnel link changes at each restart, there is no uptime guarantee, and Cloudflare quick tunnels are intended for testing and demos.
- The Mac sleeping, losing Wi-Fi or running out of disk stops the service. Docking competes with your own work for CPU.
