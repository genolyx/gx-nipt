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
    openjdk-11-jdk default-jre \
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
RUN mkdir -p /Work/NIPT/bin /Work/NIPT/data /Work/NIPT/data/bed /Work/NIPT/analysis /Work/NIPT/log /Work/NIPT/fastq

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
    pysam==0.22.0

RUN /opt/conda/envs/nipt/bin/pip install -U git+https://github.com/CenterForMedicalGeneticsGhent/WisecondorX && \
    /opt/conda/envs/nipt/bin/wisecondorx --help || echo "WisecondorX installed"

# Symlinks for convenience
RUN if [ ! -f /usr/local/bin/wisecondorx ]; then ln -s /opt/conda/envs/nipt/bin/wisecondorx /usr/local/bin/wisecondorx; fi && \
    if [ ! -f /usr/local/bin/readCounter ]; then ln -s /opt/conda/envs/nipt/bin/readCounter /usr/local/bin/readCounter; fi

# ----------------------------
# Wisecondor (Python 2) - original
# ----------------------------
RUN git clone https://github.com/VUmcCGP/wisecondor.git /opt/wisecondor && \
    sed -i '1i #!/usr/bin/env python2' /opt/wisecondor/wisecondor.py && \
    chmod +x /opt/wisecondor/wisecondor.py && \
    ln -s /opt/wisecondor/wisecondor.py /usr/local/bin/wisecondor

# ----------------------------
# Wisecondor-FF (local)
# ----------------------------
COPY bin/scripts/WisecondorFF /opt/WisecondorFF
RUN chmod +x /opt/WisecondorFF/wisecondorff && \
    sed -i 's|/opt/WisecondorFF/main.py|/opt/WisecondorFF/src/main.py|g' /opt/WisecondorFF/wisecondorff && \
    ln -s /opt/WisecondorFF/wisecondorff /usr/local/bin/wisecondor-ff

# HMMcopy reference WIG files
RUN mkdir -p /opt/conda/envs/nipt/share/hmmcopy
COPY data/wig/hg19.gc.50kb.wig  /opt/conda/envs/nipt/share/hmmcopy/
COPY data/wig/hg19.map.50kb.wig /opt/conda/envs/nipt/share/hmmcopy/
COPY data/wig/hg19.gc.10mb.wig  /opt/conda/envs/nipt/share/hmmcopy/
COPY data/wig/hg19.map.10mb.wig /opt/conda/envs/nipt/share/hmmcopy/

# AWK/R scripts
COPY bin/scripts/fetus.awk /Work/NIPT/bin/
COPY bin/scripts/mom.awk   /Work/NIPT/bin/
COPY bin/scripts/HMMcopy.R /Work/NIPT/bin/

# ----------------------------
# SeqFF setup
# ----------------------------
RUN mkdir -p /opt/nipt/bin/scripts/seqFF_R
COPY bin/scripts/seqFF_R /opt/nipt/bin/scripts/seqFF_R
RUN if [ -f /opt/nipt/bin/scripts/seqFF_R/seqff.r ]; then \
      sed -i 's|load("pd4615-sup-0008-file1.rdata")|load("/opt/nipt/bin/scripts/seqFF_R/pd4615-sup-0008-file1.rdata")|g' \
        /opt/nipt/bin/scripts/seqFF_R/seqff.r; \
    fi

# SeqFF wrapper
RUN printf '#!/bin/bash\nexport PATH=/opt/conda/envs/nipt/bin:$PATH\ncd /opt/nipt/bin/scripts/seqFF_R\n/opt/conda/envs/nipt/bin/Rscript seqff.r "$@"\n' \
    > /usr/local/bin/seqff && chmod +x /usr/local/bin/seqff

# ----------------------------
# Pipeline modules & scripts (코드 변경은 여기서부터 캐시 무효)
# ----------------------------
ARG CACHE_BREAK=1

COPY bin/scripts/modules/prizm_runner.py            /Work/NIPT/bin/
COPY bin/scripts/modules/ezd_runner.py              /Work/NIPT/bin/
COPY bin/scripts/modules/sca_detector.py            /Work/NIPT/bin/
COPY bin/scripts/modules/process_md_result.py       /Work/NIPT/bin/
COPY bin/scripts/modules/process_md_result_dev.py   /Work/NIPT/bin/
COPY bin/scripts/modules/generate_json_output.py    /Work/NIPT/bin/
COPY bin/scripts/modules/html_review_page.py        /Work/NIPT/bin/
COPY bin/scripts/modules/gender_detector.py         /Work/NIPT/bin/

COPY bin/scripts/nipt_pipeline.py /Work/NIPT/bin/
RUN chmod +x /Work/NIPT/bin/nipt_pipeline.py

# MD-only pipeline
COPY bin/scripts/md_pipeline.py /Work/NIPT/bin/scripts/
RUN chmod +x /Work/NIPT/bin/scripts/md_pipeline.py

# WC family ref generator
COPY bin/scripts/generate_wcfamily_reference.sh /Work/NIPT/bin/
RUN chmod +x /Work/NIPT/bin/generate_wcfamily_reference.sh

# for New WisecondorX build
COPY data/empty_regions.txt /Work/NIPT/data

# ----------------------------
# Env vars
# ----------------------------
ENV SRC_DIR=/Work/NIPT/bin \
    ANALYSIS_DIR=/Work/NIPT/analysis \
    LOG_DIR=/Work/NIPT/log \
    FASTQ_DIR=/Work/NIPT/fastq \
    BWA2=bwa-mem2 \
    SAMTools=samtools \
    PICARD=picard \
    qualimap=qualimap \
    WC=/opt/wisecondor \
    WCX=wisecondorx \
    WCFF=wisecondor-ff \
    HMMcopy=/opt/conda/envs/nipt \
    PYTHON2=python2.7 \
    Rscript=/opt/conda/envs/nipt/bin/Rscript \
    R_HOME=/opt/conda/envs/nipt/lib/R \
    R_LIBS_USER=/opt/conda/envs/nipt/lib/R/library \
    PATH=/Work/NIPT/bin:/opt/conda/envs/nipt/bin:$PATH

# conda env auto-activate on shell
RUN echo "source /opt/conda/etc/profile.d/conda.sh && conda activate nipt" >> /root/.bashrc

WORKDIR /Work/NIPT

# ----------------------------
# Entrypoint
# ----------------------------
ENTRYPOINT ["python3", "/Work/NIPT/bin/nipt_pipeline.py"]
CMD ["--help"]
