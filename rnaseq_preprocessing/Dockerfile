FROM nfcore/rnaseq:1.4.2
WORKDIR /install
# Update apt and install unzip
RUN apt-get update && apt-get install unzip
RUN apt-get install -y apt-utils
# Install snakemake and dependencies
RUN apt install -y build-essential
#RUN apt-get install -y gcc 
RUN apt-get install -y python3-dev
RUN  apt-get install -y python3-pip
RUN python3 -m pip install snakemake
# Install pandas library
RUN python3 -m pip install pandas==1.1.2 && python3 -m pip install gffutils==0.10.1
RUN apt-get install bedtools
#RUN apt install -y build-essential
RUN apt install zlib1g-dev
RUN wget https://github.com/alexdobin/STAR/archive/2.7.7a.tar.gz && tar -zxvf 2.7.7a.tar.gz && cd STAR-2.7.7a/source/ && make STAR

WORKDIR /work