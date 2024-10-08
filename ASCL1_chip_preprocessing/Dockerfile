FROM nfcore/chipseq:1.2.1
WORKDIR /install
# Update apt and install unzip
RUN apt-get update && apt-get install unzip
# Install bowtie2
RUN wget https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.4.1/bowtie2-2.4.1-linux-x86_64.zip && \
unzip bowtie2-2.4.1-linux-x86_64.zip && cp bowtie2-2.4.1-linux-x86_64/bowtie2 /usr/bin && cp bowtie2-2.4.1-linux-x86_64/bowtie2-align* /usr/bin && \
cp bowtie2-2.4.1-linux-x86_64/bowtie2-build* /usr/bin && cp bowtie2-2.4.1-linux-x86_64/bowtie2-inspect* /usr/bin && rm bowtie2-2.4.1-linux-x86_64.zip
# Install snakemake and dependencies
RUN apt-get install -y gcc python3-dev && pip install snakemake
# Install pandas library
RUN pip3 install pandas
# Install pysam and tidyverse
RUN pip install pysam
RUN R -e "install.packages('tidyverse', repos='http://cran.rstudio.com/')"
# Install ROSE
RUN git clone https://bitbucket.org/young_computation/rose.git
RUN conda install -c bioconda snakemake-wrapper-utils
RUN conda create -n python2 python=2.7
WORKDIR /work
COPY . /
