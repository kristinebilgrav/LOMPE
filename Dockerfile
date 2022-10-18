FROM ubuntu:18.04

RUN apt-get update && apt-get install -y python3.7 python3-pip git zlib1g-dev perl gcc make libdbi-perl wget
RUN pip3 install --upgrade pip
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-py38_4.10.3-Linux-x86_64.sh && bash Miniconda3-py38_4.10.3-Linux-x86_64.sh -p /miniconda -b

ENV PATH=/miniconda/bin:$PATH

RUN conda update -y conda \
    && rm Miniconda3-py38_4.10.3-Linux-x86_64.sh

RUN conda config --add channels defaults 
RUN conda config --add channels conda-forge 
RUN conda config --add channels bioconda 

RUN conda install -c bioconda pysam bcftools delly 

#SVDB
RUN  pip install numpy SVDB

#Nanopolish
RUN pip install -r scripts/requirements.txt --user
RUN git clone --recursive https://github.com/jts/nanopolish.git \
    && cd nanopolish \ 
    && make

