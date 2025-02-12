# Base image for SQANTI3/v5.3.0 with Ubuntu 22.04

# Using ubuntu 22.04
# Right now edlib doesn't work with python 3.12 which is the default version
# of python in Ubuntu 24.04. edlib has had no updates since April 19, 2023
# so no compatibility is expected for the time being.
# Dockerfile originally developed by @skchronicles, updated and optimized
# for version v5.3.0 by @Fabian-RY

FROM ubuntu:22.04
SHELL ["/bin/bash", "--login" ,"-c"]

# Set the versions of different softwares dependencies and SQANTI3 version
# To install
ENV SQANTI3_VERSION="5.3.6"
ENV DESALT_VERSION="1.5.6"
ENV NAMFINDER_VERSION="0.1.3"

LABEL maintainer="aarzalluz" \
   base_image="ubuntu:22.04" \
   version="v0.9.0"   \
   software="sqanti3/v${SQANTI3_VERSION}" \
   about.summary="SQANTI3: Tool for the Quality Control of Long-Read Defined Transcriptomes" \
   about.home="https://github.com/ConesaLab/SQANTI3" \
   about.documentation="https://github.com/ConesaLab/SQANTI3/wiki/" \
   about.tags="Transcriptomics"

# Default Timeset inside the container: UTC

############### INIT ################
# Create Container filesystem specific 
# working directory and opt directories
# to avoid collisions with the host's
# filesystem, i.e. /opt and /data
RUN mkdir -p /opt2 && mkdir -p /data2
WORKDIR /opt2


############### SETUP ################
# This section installs system packages 
# required for your project. If you need 
# extra system packages add them here.
RUN apt-get update \
   && apt-get -y upgrade \
   && DEBIAN_FRONTEND=noninteractive apt-get install -y \
       # bedtools/2.30.0
       bedtools \
       build-essential \
       cmake \
       cpanminus \
       curl \
       gawk \
       # gffread/0.12.7
       gffread \
       # gmap/2021-12-17
       gmap \
       gzip \
       # kallisto/0.46.2
       kallisto \
       libcurl4-openssl-dev \
       libssl-dev \
       libxml2-dev \
       locales \
       # minimap2/2.24
       minimap2 \
       # perl/5.34.0-3
       perl \
       pkg-config \
       # STAR/2.7.10a
       rna-star \
       # samtools/1.13-4
       samtools \
       # seqtk/1.3-2
       seqtk \
       wget \
       zlib1g-dev \
   && apt-get clean && apt-get purge \
   && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# Set the locale
RUN localedef -i en_US -f UTF-8 en_US.UTF-8
# Perl fix issue
RUN cpanm FindBin Term::ReadLine

## Installing miniconda inside the container

############### MANUAL ################
# Install tools from src manually,
# Installs deSALT from GitHub: (current version, 1.5.6, defined in ENV variable)
# https://github.com/ydLiu-HIT/deSALT/releases/tag/v1.5.6
# This tool was created using an older
# version of GCC that allowed multiple
# definitions of global variables.
# We are using GCC/10, which does not
# allow multiple definitions. Adding
# -Wl,--allow-multiple-definition
# to the linker to fix this issue.
RUN mkdir -p /opt2/desalt/${DESALT_VERSION}/ \
   && wget https://github.com/ydLiu-HIT/deSALT/archive/refs/tags/v${DESALT_VERSION}.tar.gz -O /opt2/desalt/${DESALT_VERSION}/v${DESALT_VERSION}.tar.gz \
   && tar -zvxf /opt2/desalt/${DESALT_VERSION}/v${DESALT_VERSION}.tar.gz -C /opt2/desalt/${DESALT_VERSION}/ \
   && rm -f /opt2/desalt/${DESALT_VERSION}/v${DESALT_VERSION}.tar.gz \
   && cd /opt2/desalt/${DESALT_VERSION}/deSALT-${DESALT_VERSION}/src/deBGA-master/ \
   && make CFLAGS="-g -Wall -O2 -Wl,--allow-multiple-definition" \
   && cd .. \
   && make CFLAGS="-g -Wall -O3 -Wc++-compat -Wno-unused-variable -Wno-unused-but-set-variable -Wno-unused-function -Wl,--allow-multiple-definition"

ENV PATH="${PATH}:/opt2/desalt/${DESALT_VERSION}/deSALT-${DESALT_VERSION}/src"
WORKDIR /opt2

# Installs namfinder, requirement of
# ultra-bioinformatics tool from pypi.
RUN mkdir -p /opt2/namfinder/${NAMFINDER_VERSION}/ \
   && wget https://github.com/ksahlin/namfinder/archive/refs/tags/v${NAMFINDER_VERSION}.tar.gz -O /opt2/namfinder/${NAMFINDER_VERSION}/v${NAMFINDER_VERSION}.tar.gz \
   && tar -zvxf /opt2/namfinder/${NAMFINDER_VERSION}/v${NAMFINDER_VERSION}.tar.gz -C /opt2/namfinder/${NAMFINDER_VERSION}/ \
   && rm -f /opt2/namfinder/${NAMFINDER_VERSION}/v${NAMFINDER_VERSION}.tar.gz \
   && cd /opt2/namfinder/${NAMFINDER_VERSION}/namfinder-${NAMFINDER_VERSION}/ \
   # Build to be compatiable with most
   # Intel x86 CPUs, should work with
   # old hardware, i.e. sandybridge
   && cmake -B build -DCMAKE_C_FLAGS="-msse4.2" -DCMAKE_CXX_FLAGS="-msse4.2" \
   && make -j -C build

ENV PATH="${PATH}:/opt2/namfinder/${NAMFINDER_VERSION}/namfinder-${NAMFINDER_VERSION}/build"
WORKDIR /opt2

########### SQANTI3 (currentily v${NAMFINDER_VERSION}) ############
# Installs SQANTI3 with the version defined in the ENV variable (currently 5.3.0)
# dependenciesand requirements have already been
# satisfied, for more info see:
# https://github.com/ConesaLab/SQANTI3
RUN mkdir -p /opt2/sqanti3/${SQANTI3_VERSION}/ \
   && wget https://github.com/ConesaLab/SQANTI3/archive/refs/tags/v${SQANTI3_VERSION}.tar.gz -O /opt2/sqanti3/${SQANTI3_VERSION}/v${SQANTI3_VERSION}.tar.gz \
   && tar -zvxf /opt2/sqanti3/${SQANTI3_VERSION}/v${SQANTI3_VERSION}.tar.gz -C /opt2/sqanti3/${SQANTI3_VERSION}/ \
   && rm -f /opt2/sqanti3/${SQANTI3_VERSION}/v${SQANTI3_VERSION}.tar.gz \
   # Removing exec bit for non-exec files
   && chmod -x \
       /opt2/sqanti3/${SQANTI3_VERSION}/SQANTI3-${SQANTI3_VERSION}/LICENSE \
       /opt2/sqanti3/${SQANTI3_VERSION}/SQANTI3-${SQANTI3_VERSION}/.gitignore \
       /opt2/sqanti3/${SQANTI3_VERSION}/SQANTI3-${SQANTI3_VERSION}/*.md \
       /opt2/sqanti3/${SQANTI3_VERSION}/SQANTI3-${SQANTI3_VERSION}/*.yml \
   # Patch: adding absolute PATH to howToUse.png
   # that gets embedded in the report. When running
   # sqanti_qc.py within docker/singularity container,
   # it fails at the report generation step because 
   # pandoc cannot find the png file (due to relative
   # path). Converting relative path in Rmd files to
   # an absolute path to avoid this issue altogether.
   && sed -i \
       's@src="howToUse.png"@src="/opt2/sqanti3/${SQANTI3_VERSION}/SQANTI3-${SQANTI3_VERSION}/src/utilities/report_qc/howToUse.png"@g' \
       /opt2/sqanti3/${SQANTI3_VERSION}/SQANTI3-${SQANTI3_VERSION}/src/utilities/report_qc/SQANTI3_report.Rmd \
       /opt2/sqanti3/${SQANTI3_VERSION}/SQANTI3-${SQANTI3_VERSION}/src/utilities/report_pigeon/pigeon_report.Rmd

ENV PATH="${PATH}:/opt2/sqanti3/${SQANTI3_VERSION}/SQANTI3-${SQANTI3_VERSION}:/opt2/sqanti3/${SQANTI3_VERSION}/SQANTI3-${SQANTI3_VERSION}/src/utilities"
WORKDIR /opt2/sqanti3/${SQANTI3_VERSION}/SQANTI3-${SQANTI3_VERSION}

RUN mkdir -p /conda/miniconda3 && \
    wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /conda/miniconda3/miniconda.sh && \
    bash /conda/miniconda3/miniconda.sh -b -u -p /conda/miniconda3 && \
    rm /conda/miniconda3/miniconda.sh && \
    /conda/miniconda3/bin/conda env create -f SQANTI3.conda_env.yml && \
    /conda/miniconda3/bin/conda clean -a


################ POST #################
# Add Dockerfile and export environment 
# variables and update permissions
ADD Dockerfile /opt2/sqanti3_${SQANTI3_VERSION}.dockerfile
ENV PATH="/opt2:$PATH"
# Hide deprecation warnings from sqanit
ENV PYTHONWARNINGS="ignore::DeprecationWarning"
WORKDIR /data2

ENV PATH="${PATH}:/conda/miniconda3/bin/"

RUN chmod -R a+rX /opt2 && \
    ln -s /opt2/sqanti3/${SQANTI3_VERSION}/SQANTI3-${SQANTI3_VERSION}/sqanti3_qc.py sqanti3_qc.py && \
    ln -s /opt2/sqanti3/${SQANTI3_VERSION}/SQANTI3-${SQANTI3_VERSION}/sqanti3_filter.py sqanti3_filter.py && \
    ln -s /opt2/sqanti3/${SQANTI3_VERSION}/SQANTI3-${SQANTI3_VERSION}/sqanti3_rescue.py sqanti3_rescue.py && \
    apt-get remove -y build-essential cmake && apt-get autoremove -y && apt-get clean -y && \
    apt-get clean autoclean -y

# Creating an user for non root executions inside the container...
RUN groupadd user && useradd -r -g user user && chown user:user .

# ... and efectively using that user
USER user

ENTRYPOINT ["conda", "run", "--no-capture-output" ,"-n","SQANTI3.env"]
CMD ["/bin/bash"]
