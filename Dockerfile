FROM continuumio/miniconda3

# Environment
ENV DEBIAN_FRONTEND=noninteractive \
    LANG=en_US.UTF-8 \
    LC_ALL=C.UTF-8 \
    LANGUAGE=en_US.UTF-8

RUN apt-get update --fix-missing && \
    apt-get install \
        --no-install-recommends \
        --no-install-suggests \
        --fix-broken \
        --quiet \
        --assume-yes \
        --no-show-upgraded \
            wget \
            gcc \
            gfortran \
            git \
            g++ \
            libcurl4-openssl-dev \
            libigraph-dev \
            libssl-dev \
            libxml2-dev \
            minimap2 \
            pandoc \
            procps \
            locales \
            localepurge && \
    locale-gen en_US.UTF-8 && \
    localepurge && \
    apt-get clean && \
    rm -rf /tmp/downloaded_packages/* && \
    apt-get purge && \
    apt-get autoclean

WORKDIR /app

# Create the environment:
COPY SQANTI3.conda_env.yml .
RUN conda env create -f SQANTI3.conda_env.yml

# Make RUN commands use the new environment:
RUN echo "conda activate SQANTI3.env" >> ~/.bashrc
SHELL ["/bin/bash", "--login", "-c"]

RUN git clone https://github.com/Magdoll/cDNA_Cupcake.git && \
    cd cDNA_Cupcake && \
    python setup.py build && \
    python setup.py install

# Put conda in path so we can use conda activate
ENV PYTHONPATH=/app/cDNA_Cupcake/sequence:$PYTHONPATH

# Copy rest of source
COPY . .

# The code to run when container is started:
RUN chmod +x entrypoint.sh
ENTRYPOINT ["./entrypoint.sh"]
