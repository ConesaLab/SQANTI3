# translated from https://github.com/mestia/SQANTI3/blob/master/SQANTI3_singularity_recipe.def
FROM debian:buster-slim

ENV DEBIAN_FRONTEND=noninteractive

RUN apt update && \
    apt install \
        -y \
        --no-install-recommends \
        --no-install-suggests \
            localepurge \
            build-essential \
            git \
            ca-certificates \
            wget \
            curl \
            libkrb5-3 \
            libk5crypto3 && \
    # the deps for gtfToGenePred, see ldd gtfToGenePred
    locale-gen en_US.UTF-8 && \
    localepurge && \
    apt-get clean && \
    apt-get autoclean

ADD sqanti3_env.yml /opt/sqanti3_env.yml
ADD https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh /opt/Miniconda3.sh
RUN bash /opt/Miniconda3.sh -b -p /opt/miniconda3 && \
    # git clone https://github.com/milescsmith/SQANTI3.git /opt/SQANTI3 && \
    eval "$(/opt/miniconda3/bin/conda shell.bash hook)" && \
    # it seems this yaml need some fresh packages
    conda env create -f /opt/sqanti3_env.yml && \
    conda init bash && \
    conda clean -a

ENV LC_ALL=C
ENV PATH=/opt/miniconda3/envs/sqanti3/bin:/opt/miniconda3/envs/sqanti3/utilities:/SQANTI3:$PATH
ENV PYTHONPATH=/opt/miniconda3/envs/sqanti3/lib/python3.7/site-packages/

# ADD http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/gtfToGenePred /opt/gtfToGenePred/gtfToGenePred
# ENV PATH=/opt/gtfToGenePred:$PATH

RUN /opt/miniconda3/envs/sqanti3/bin/pip install git+https://github.com/milescsmith/cDNA_Cupcake.git@prune
RUN /opt/miniconda3/envs/sqanti3/bin/pip install git+https://github.com/milescsmith/pygmst.git

COPY . /opt/sqanti3
RUN /opt/miniconda3/envs/sqanti3/bin/pip install -e /opt/sqanti3
# cleanup
RUN apt clean && \
    apt autoremove -y && \
    # rm -rf /include
    rm -rf /opt/Miniconda3-latest-Linux-x86_64.sh && \
    rm -rf /opt/cDNA_Cupcake/.git && \
    rm -rf /opt/SQANTI3/.git && \
    rm -fr /opt/sqanti3_env.yml

ENTRYPOINT [ "sqanti3_qc" ]
