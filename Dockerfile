# pscreen — hosted (public web) image. Linux x86-64, CPU only. Built for Hugging Face Docker Spaces, works on any Docker host.
# NOT yet built or run by the author on Linux (the development machine is a Mac without Docker): the first build is the first test.
FROM python:3.10-slim-bookworm

ENV DEBIAN_FRONTEND=noninteractive PIP_NO_CACHE_DIR=1 PYTHONUNBUFFERED=1 \
    PSCREEN_HOSTED=1 PSCREEN_DATA=/data PSCREEN_SLOTS=1 PSCREEN_RETENTION_H=48 \
    P2RANK_HOME=/app/third_party/p2rank_2.5 HOME=/home/user
RUN apt-get update && apt-get install -y --no-install-recommends openjdk-17-jre-headless libgomp1 libxrender1 libxext6 libsm6 curl git ca-certificates unzip \
    && rm -rf /var/lib/apt/lists/* && useradd -m -u 1000 user && mkdir -p /app /data && chown -R user /app /data
USER user
WORKDIR /app

# main environment (UI, RDKit, analysis)
RUN pip install --user rdkit pandas numpy scipy streamlit biopython openpyxl
ENV PATH=/home/user/.local/bin:$PATH

# P2Rank 2.5 (MIT)
RUN mkdir -p third_party && curl -fsSL https://github.com/rdk/p2rank/releases/download/2.5/p2rank_2.5.tar.gz | tar xz -C third_party

# GNINA 1.3.3 (Apache-2.0) + the CUDA runtime libraries the static binary needs in order to load (used with --no_gpu)
RUN mkdir -p bin infra/cudalibs && curl -fsSL -o bin/gnina.cuda12.8.static https://github.com/gnina/gnina/releases/download/v1.3.3/gnina.cuda12.8.static && chmod +x bin/gnina.cuda12.8.static \
    && pip install --target /tmp/cu nvidia-cuda-runtime-cu12 nvidia-cudnn-cu12 nvidia-cublas-cu12 nvidia-cusparse-cu12 nvidia-cufft-cu12 nvidia-cusolver-cu12 nvidia-nvjitlink-cu12 \
    && find /tmp/cu -name "*.so*" -exec cp -P {} infra/cudalibs/ \; && rm -rf /tmp/cu

# PocketMiner (MIT per its README) in its own environment: CPU TensorFlow 2.10 + the pinned versions from docs/SETUP.md
RUN git clone https://github.com/bowman-lab/ae-pocketminer third_party/ae-pocketminer && git -C third_party/ae-pocketminer checkout 3328ca0 \
    && python -m venv envs/pocketminer && envs/pocketminer/bin/pip install tensorflow-cpu==2.10.0 numpy==1.26.4 mdtraj==1.10.0 protobuf==3.19.6 pyyaml scipy tqdm scikit-learn pandas h5py

# Lacuna 1.2.0 (MIT) and ADMET-AI 2.0.1 (MIT), each isolated
RUN python -m venv envs/lacuna && envs/lacuna/bin/pip install lacuna-pockets \
    && python -m venv envs/admet && envs/admet/bin/pip install admet-ai

COPY --chown=user . /app
EXPOSE 7860
CMD ["sh", "-c", "streamlit run app.py --server.port ${PORT:-7860} --server.address 0.0.0.0 --server.headless true --server.enableXsrfProtection true --server.maxUploadSize 20 --browser.gatherUsageStats false"]
