FROM ubuntu:16.04

# install all build-essentials
RUN apt-get update \
	&& apt-get install -y --no-install-recommends \
        bc \
        bzip2 \
        ca-certificates \
        curl \
        g++ \
        gcc \
        gfortran \
        git \
        less \
        libbz2-dev \
        libcurl4-openssl-dev \
        libgsl-dev \
        libgsl2 \
        liblzma-dev \
        libmysqlclient-dev \
        libncurses5-dev \
        libpcre3-dev \
        libreadline-dev \
        libssl-dev \
        libx11-dev \
        make \
        python-dev \
        texlive-latex-base texlive-latex-extra \
        unzip \
        xorg-dev \
        zlib1g-dev \
    && rm -rf /var/lib/apt/lists/* 

# install R 3.3.2
RUN cd /opt/ \
    && curl -fsSL https://repo.continuum.io/miniconda/Miniconda2-latest-Linux-x86_64.sh -o /opt/Miniconda2-latest-Linux-x86_64.sh \
    && bash Miniconda2-latest-Linux-x86_64.sh -b -p /opt/miniconda \
    && rm Miniconda2-latest-Linux-x86_64.sh \
    && ln -s /opt/miniconda/bin/conda /usr/bin/ \
    && conda config --add channels defaults \
    && conda config --add channels conda-forge \
    && conda config --add channels bioconda \
    && conda install r=3.3.2 \
    && echo 'install.packages("ggplot2",repos="http://cran.us.r-project.org")' > /opt/packages.R \
    && echo 'install.packages("cowplot",repos="http://cran.us.r-project.org")' >> /opt/packages.R \
    && echo 'install.packages("jsonlite",repos="http://cran.us.r-project.org")' >> /opt/packages.R \
    && /opt/miniconda/bin/Rscript /opt/packages.R \
    && rm /opt/packages.R \
    && mkdir -p /usr/local/lib/R/site-library

# install fastq-dump.2.8.2
RUN curl -fsSL http://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-ubuntu64.tar.gz -o /opt/sratoolkit.current-ubuntu64.tar.gz \
	&& tar -xvzf /opt/sratoolkit.current-ubuntu64.tar.gz -C /opt/ \
	&& rm /opt/sratoolkit.current-ubuntu64.tar.gz

# install fastqc
RUN conda install -c bioconda fastqc

# install preseq 2.0.0
RUN curl -fsSL http://smithlabresearch.org/downloads/preseq_linux_v2.0.tar.bz2 -o /opt/preseq_linux_v2.0.tar.bz2 \
    && tar xvjf /opt/preseq_linux_v2.0.tar.bz2 -C /opt/ \
    && ln -s /opt/preseq_v2.0/preseq /usr/local/bin/preseq \
    && ln -s /opt/preseq_v2.0/bam2mr /usr/local/bin/bam2mr \
    && rm /opt/preseq_linux_v2.0.tar.bz2 \
    # make sure that libgsl.so.0 exists beause PreSeq links to that
    && ln -s /usr/lib/x86_64-linux-gnu/libgsl.so /lib/x86_64-linux-gnu/libgsl.so.0

# install cutadapt 1.14
RUN curl -fsSL https://bootstrap.pypa.io/get-pip.py -o /opt/get-pip.py \
    && python /opt/get-pip.py \
    && rm /opt/get-pip.py \
    && pip install cutadapt

# install bwa 0.7.16a
RUN curl -fsSL https://sourceforge.net/projects/bio-bwa/files/bwa-0.7.16a.tar.bz2 -o /opt/bwa-0.7.16a.tar.bz2 \
    && tar xvjf /opt/bwa-0.7.16a.tar.bz2 -C /opt/ \
	&& cd /opt/bwa-0.7.16a/ \
	&& make \
	&& rm /opt/bwa-0.7.16a.tar.bz2

# install macs2 v2.1.1.20160309
RUN pip install numpy scipy \
	&& pip install MACS2

# install samtools v1.2
RUN curl -fsSL https://github.com/samtools/samtools/releases/download/1.3.1/samtools-1.3.1.tar.bz2 -o /opt/samtools-1.3.1.tar.bz2 \
    && tar xvjf /opt/samtools-1.3.1.tar.bz2 -C /opt/ \
    && cd /opt/samtools-1.3.1 \
    && make \
    && make install \
    && rm /opt/samtools-1.3.1.tar.bz2

# install bedtools 2.25.0
RUN curl -fsSL https://github.com/arq5x/bedtools2/releases/download/v2.25.0/bedtools-2.25.0.tar.gz -o /opt/bedtools-2.25.0.tar.gz \
	&& tar -zxvf /opt/bedtools-2.25.0.tar.gz -C /opt/ \
	&& cd /opt/bedtools2 \
	&& make \
	&& rm /opt/bedtools-2.25.0.tar.gz

# install methylQA 0.2.0
RUN cd /opt/ \
	&& git clone git://github.com/lidaof/methylQA.git \
	&& cd /opt/methylQA/ \
	&& make \
	&& cd /

# install bedGraphToBigWig
RUN conda install ucsc-bedgraphtobigwig
    
ENV PATH /opt/sratoolkit.2.8.2-1-ubuntu64/bin:/opt/bwa-0.7.16a/:/opt/bedtools2/bin/:/opt/methylQA/:/opt/kentUtils/bin/:/opt/miniconda/bin/:$PATH
