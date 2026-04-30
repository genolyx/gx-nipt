# syntax=docker/dockerfile:1.7
# Base
FROM ubuntu:22.04

# ----------------------------
# Locale / timezone / apt setup
# ----------------------------
ENV DEBIAN_FRONTEND=noninteractive
ENV TZ=Asia/Seoul

RUN apt-get update && \
    apt-get install -y tzdata && \
    ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

RUN apt-get update && \
    apt-get install -y --no-install-recommends software-properties-common && \
    add-apt-repository universe && \
    apt-get update && apt-get clean && \
    rm -rf /var/lib/apt/lists/* /var/cache/apt/archives/*

# ----------------------------
# OS packages (bio + build deps)
# ----------------------------
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    python3 python3-pip \
    python2 python2-dev \
    curl \
    openjdk-17-jre-headless \
    git wget unzip \
    zlib1g-dev libbz2-dev liblzma-dev libncurses5-dev \
    libcurl4-openssl-dev libssl-dev libxml2-dev \
    libncurses-dev libreadline-dev \
    gfortran gcc make \
    pkg-config \
    libicu-dev libicu70 \
    bedtools \
    bwa samtools \
    fastqc seqtk \
    vim zip \
    && apt-get clean && \
    rm -rf /var/lib/apt/lists/* /var/cache/apt/archives/*

# Python 2 pip + py2 deps (Wisecondor original)
RUN curl https://bootstrap.pypa.io/pip/2.7/get-pip.py | python2.7 && \
    python2.7 -m pip install numpy==1.14.6 matplotlib pysam scipy futures subprocess32 scikit-learn==0.20.4

# Host python3 pip upgrade (별도)
RUN python3 -m pip install --upgrade pip

# Workdirs
# /Work/NIPT/* is kept only for legacy tools (picard jar, conda env); all
# pipeline-owned scripts now live under /opt/gx-nipt/bin/scripts/ and are
# called directly by Nextflow modules.
RUN mkdir -p /Work/NIPT/bin /Work/NIPT/analysis /Work/NIPT/log /Work/NIPT/fastq \
             /opt/gx-nipt/bin/scripts

# ----------------------------
# External tool: bwa-mem2
# ----------------------------
RUN cd /tmp && \
    wget https://github.com/bwa-mem2/bwa-mem2/releases/download/v2.2.1/bwa-mem2-2.2.1_x64-linux.tar.bz2 && \
    tar -xvjf bwa-mem2-2.2.1_x64-linux.tar.bz2 && \
    cp bwa-mem2-2.2.1_x64-linux/bwa-mem2* /usr/local/bin/ && \
    ln -sf /usr/local/bin/bwa-mem2.avx2 /usr/local/bin/bwa-mem2 && \
    rm -rf /tmp/*

# ----------------------------
# Picard (jar + wrapper)
# ----------------------------
RUN mkdir -p /Work/NIPT/bin/picard && \
    wget https://github.com/broadinstitute/picard/releases/download/2.27.5/picard.jar -O /Work/NIPT/bin/picard/picard.jar && \
    echo '#!/bin/bash\njava "$@" -jar /Work/NIPT/bin/picard/picard.jar' > /usr/local/bin/picard && \
    chmod +x /usr/local/bin/picard

# ----------------------------
# Miniconda (single) + pin .condarc (conda-forge/bioconda only)
# ----------------------------
ENV PATH=/opt/conda/bin:$PATH

RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /tmp/miniconda.sh && \
    bash /tmp/miniconda.sh -b -p /opt/conda && rm /tmp/miniconda.sh && \
    conda clean -afy

# .condarc 고정: defaults 완전 제거, strict priority
RUN printf "%s\n" \
"channels:" \
"  - conda-forge" \
"  - bioconda" \
"channel_priority: strict" \
"default_channels: []" \
"show_channel_urls: true" \
> /root/.condarc && cp /root/.condarc /opt/conda/.condarc || true

# libmamba solver 설치 (conda-forge만 사용)
RUN conda install -n base -y conda-libmamba-solver --override-channels -c conda-forge && \
    conda config --set solver libmamba

# ----------------------------
# Create env (Py 3.11.9) - override channels + repodata.json
# ----------------------------
RUN --mount=type=cache,target=/opt/conda/pkgs \
    conda create -y -n nipt python=3.11.9 \
      --override-channels -c conda-forge \
      --repodata-fn=repodata.json

# ----------------------------
# Bio/R stack install (one-shot solve)
# ----------------------------
RUN --mount=type=cache,target=/opt/conda/pkgs \
    bash -lc '\
      source /opt/conda/etc/profile.d/conda.sh && conda activate nipt && \
      conda install -y \
        --override-channels -c conda-forge -c bioconda --strict-channel-priority \
        --repodata-fn=repodata.json \
        hmmcopy \
        qualimap \
        r-base=4.4 \
        r-essentials \
        bioconductor-rsamtools \
        bioconductor-genomicalignments \
        bioconductor-genomicranges \
        bioconductor-rtracklayer \
        bioconductor-s4vectors \
        bioconductor-iranges \
        bioconductor-hmmcopy \
        bioconductor-dnacopy \
        r-argparse r-optparse r-r6 r-jsonlite r-data.table r-ggplot2 r-plyr r-mass r-e1071 \
        r-stringi r-stringr \
      && conda clean -afy'

# R stringr quick check
RUN /opt/conda/envs/nipt/bin/R --slave -e "\
  library(stringr); \
  result <- str_replace('test', 'test', 'success'); \
  if (result != 'success') stop('StringR test failed'); \
  message('✓ StringR verified successfully');"

# seqFF 보완 (CRAN)
RUN /opt/conda/envs/nipt/bin/R --slave -e "\
  install.packages(c('tibble','dplyr'), repos='https://cloud.r-project.org', dependencies=TRUE, force=TRUE)"

# ----------------------------
# WisecondorX deps + install (pip in env)
# ----------------------------
RUN /opt/conda/envs/nipt/bin/pip install \
    numpy==1.26.4 \
    scipy==1.13.0 \
    scikit-learn==1.4.2 \
    plotly \
    pandas==2.2.2 \
    matplotlib==3.8.4 \
    pysam==0.22.0 \
    typer \
    loguru

RUN /opt/conda/envs/nipt/bin/pip install -U git+https://github.com/CenterForMedicalGeneticsGhent/WisecondorX && \
    /opt/conda/envs/nipt/bin/wisecondorx --help || echo "WisecondorX installed"

# CBS.R is the R script called by WisecondorX for circular binary segmentation.
# The pip-installed wheel omits the include/ package data directory, so we copy
# it manually to the location WisecondorX expects (os.path.dirname(__file__)/include/).
COPY bin/wcx_include/CBS.R /tmp/CBS.R
RUN WCXDIR=$(/opt/conda/envs/nipt/bin/python -c \
      "import wisecondorx, os; print(os.path.dirname(wisecondorx.__file__))") && \
    mkdir -p "${WCXDIR}/include" && \
    cp /tmp/CBS.R "${WCXDIR}/include/CBS.R" && \
    echo "CBS.R installed at ${WCXDIR}/include/CBS.R"

# Symlinks for convenience
RUN if [ ! -f /usr/local/bin/wisecondorx ]; then ln -s /opt/conda/envs/nipt/bin/wisecondorx /usr/local/bin/wisecondorx; fi && \
    if [ ! -f /usr/local/bin/readCounter ]; then ln -s /opt/conda/envs/nipt/bin/readCounter /usr/local/bin/readCounter; fi

# ----------------------------
# Nextflow (launcher container runs Nextflow itself; task containers are
# spawned on the host daemon via the bind-mounted docker socket — a.k.a.
# Docker-outside-of-Docker, same pattern as gx-exome).
# ----------------------------
RUN curl -s https://get.nextflow.io | bash && \
    mv nextflow /usr/local/bin/ && \
    chmod +x /usr/local/bin/nextflow && \
    nextflow -version

# ----------------------------
# Wisecondor (Python 2) - original
# ----------------------------
RUN git clone https://github.com/VUmcCGP/wisecondor.git /opt/wisecondor && \
    sed -i '1i #!/usr/bin/env python2' /opt/wisecondor/wisecondor.py && \
    chmod +x /opt/wisecondor/wisecondor.py && \
    ln -s /opt/wisecondor/wisecondor.py /usr/local/bin/wisecondor

# ----------------------------
# Pipeline scripts (Nextflow-called)
#
# Reference genomes / HMMcopy WIGs / lab WC(X) NPZs / ML models are NOT baked
# into the image. They are bind-mounted from the host via nextflow.config's
# docker.runOptions (-v ${params.ref_dir}:${params.ref_dir}:ro).
# See docs/reference_generation.md for the expected host layout.
#
# All scripts are placed under /opt/gx-nipt/bin/scripts/ — the single path
# that every Nextflow module imports.
# ----------------------------
ARG CACHE_BREAK=15

COPY bin/scripts /opt/gx-nipt/bin/scripts
# Old Wisecondor (Python 2) – used by RUN_WC process for WC-family analysis
COPY bin/wisecondor /opt/wisecondor
RUN chmod +x /opt/wisecondor/wisecondor.py
RUN find /opt/gx-nipt/bin/scripts -type f \( -name '*.sh' -o -name '*.py' -o -name '*.awk' -o -name '*.R' \) -print0 \
    | xargs -0 chmod +x

# ----------------------------
# Env vars
# ----------------------------
ENV GXNIPT_SCRIPTS=/opt/gx-nipt/bin/scripts \
    BWA2=bwa-mem2 \
    SAMTools=samtools \
    PICARD=picard \
    qualimap=qualimap \
    WC=/opt/wisecondor \
    WCX=wisecondorx \
    HMMcopy=/opt/conda/envs/nipt \
    PYTHON2=python2.7 \
    Rscript=/opt/conda/envs/nipt/bin/Rscript \
    R_HOME=/opt/conda/envs/nipt/lib/R \
    R_LIBS_USER=/opt/conda/envs/nipt/lib/R/library \
    PATH=/opt/gx-nipt/bin/scripts:/opt/conda/envs/nipt/bin:$PATH

# conda env auto-activate on shell
RUN echo "source /opt/conda/etc/profile.d/conda.sh && conda activate nipt" >> /root/.bashrc

WORKDIR /opt/gx-nipt

# ----------------------------
# No default ENTRYPOINT.
# Nextflow invokes specific tools inside each process; see nextflow.config's
# docker.runOptions which also includes --entrypoint "" defensively.
# For interactive use: `docker run --rm -it gx-nipt:latest bash`.
# ----------------------------
CMD ["/bin/bash"]
