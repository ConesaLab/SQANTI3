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

ENV PATH=/opt/miniconda3/envs/sqanti3/bin:/opt/miniconda3/envs/sqanti3/utilities:/SQANTI3:$PATH
ENV PYTHONPATH=/opt/miniconda3/envs/sqanti3/lib/python3.7/site-packages/

RUN git clone https://github.com/Magdoll/cDNA_Cupcake.git /opt/cDNA_Cupcake && \
    cd /opt/cDNA_Cupcake && python setup.py build && python setup.py install

ADD . /opt/SQANTI3/

ADD http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/gtfToGenePred /opt/SQANTI3/utilities/gtfToGenePred
RUN chmod +x /opt/SQANTI3/utilities/gtfToGenePred

# cleanup
RUN apt remove -y \
        build-essential \
        git \
        wget \
        curl && \
    apt autoremove -y && \
    # rm -rf /include
    rm -rf /opt/Miniconda3-latest-Linux-x86_64.sh && \
    rm -rf /opt/cDNA_Cupcake/.git && \
    rm -rf /opt/SQANTI3/.git && \
    rm -fr /opt/sqanti3_env.yml

RUN chmod +x /opt/SQANTI3/sqanti3_qc.py

ENV LC_ALL=C
ENV PATH=/opt/miniconda3/envs/SQANTI3.env/bin:/opt/miniconda3/envs/SQANTI3.env/utilities:/opt/SQANTI3:$PATH
ENV PYTHONPATH=/opt/miniconda3/envs/SQANTI3.env/lib/python3.7/site-packages:/opt/cDNA_Cupcake:/opt/cDNA_Cupcake/sequence/

ENTRYPOINT [ "/opt/SQANTI3/sqanti3_qc.py" ]