# Setup (macOS arm64, no GPU)

Main env (RDKit, pandas, Streamlit, py3Dmol, Biopython, scipy): `pip install rdkit pandas streamlit py3Dmol biopython scipy`.

## GNINA (no macOS/ARM build exists)
GNINA 1.3.3 ships Linux x86-64 binaries only. We run it in a Lima VM with Rosetta:
1. `brew install lima`; `limactl start --name=gnina infra/gnina-vm.yaml`
2. In the VM: `dpkg --add-architecture amd64`, restrict existing apt sources to arm64, add the amd64 archive, install `libc6:amd64 libstdc++6:amd64 libgomp1:amd64 zlib1g:amd64 libbz2-1.0:amd64`.
3. `gh release download v1.3.3 -R gnina/gnina -p gnina.cuda12.8.static -D bin` (2 GB).
4. Put NVIDIA's x86-64 CUDA 12 runtime libs (cudnn 9, cudart, cublas, cublasLt, cusparse, cufft, cusolver, nvjitlink; from the `nvidia-*-cu12` pip wheels) in `infra/cudalibs/`. The binary needs them to load even with `--no_gpu`.
5. `infra/gnina.sh --version`. All paths passed to GNINA must be under `~/proteinscreen` (the VM mount).

## Pocket tools (each isolated; the app reports a clear error if one is missing)
- **P2Rank 2.5** (MIT): unpack a release, set `P2RANK_HOME` (or place in `third_party/p2rank_2.5`). Needs Java 11+.
- **PocketMiner** (MIT): `git clone https://github.com/bowman-lab/ae-pocketminer third_party/ae-pocketminer` (bundles the original PocketMiner checkpoint, tested at commit 3328ca0);
  `conda create -p envs/pocketminer python=3.10; envs/pocketminer/bin/pip install tensorflow-macos==2.10.0 numpy==1.26.4 mdtraj==1.10.0 protobuf==3.19.6 pyyaml scipy tqdm scikit-learn pandas h5py`
- **Lacuna 1.2.0** (MIT): `python3 -m venv envs/lacuna; envs/lacuna/bin/pip install lacuna-pockets`

- **PyMOL 3.1 (open source)** for ray-traced renders and .pse/.pml exports: `conda create -p envs/pymol -c conda-forge pymol-open-source`
- **ADMET-AI 2.0.1** (optional property predictions, MIT; checkpoints bundled in the package): `python3 -m venv envs/admet; envs/admet/bin/pip install admet-ai`
- **Mol\*** viewer: loaded from the jsDelivr CDN (needs internet only when that viewer is selected). 3Dmol and PyMOL work offline.

## Run
`python3 -m pytest -q`; `streamlit run app.py` (theme and localhost binding come from `.streamlit/config.toml`).
CLI pipeline: `python -m pscreen.pipeline runs/<name>` (reads `run_config.json`; resumable).
