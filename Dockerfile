FROM mambaorg/micromamba:jammy-cuda-12.1.0

USER root

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y \
    wget \
    && rm -rf /var/lib/apt/lists/*

RUN wget -O /tmp/dorado.tar.gz https://cdn.oxfordnanoportal.com/software/analysis/dorado-0.7.2-linux-x64.tar.gz && \
    tar -zxvf /tmp/dorado.tar.gz &&  \
    rm /tmp/dorado.tar.gz

ENV PATH="/tmp/dorado-0.7.2-linux-x64/bin:$PATH"

USER $MAMBA_USER

WORKDIR /app

COPY --chown=$MAMBA_USER:$MAMBA_USER . /app

RUN micromamba install -y -n base -f environments/docker.yml && \
    micromamba clean --all --yes

ARG MAMBA_DOCKERFILE_ACTIVATE=1  # activate environment
RUN pip install -e .

CMD ["savanna"]

