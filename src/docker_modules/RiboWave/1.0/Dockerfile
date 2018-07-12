FROM ubuntu:16.04
MAINTAINER Emmanuel Labaronne

ENV BEDTOOLS_VERSION=2.25.0*
ENV PACKAGES bedtools=${BEDTOOLS_VERSION}* \
	r-base \
	r-base-dev \
	git

RUN apt-get update &&\ 
	apt-get install -y --no-install-recommends ${PACKAGES} && \
	apt-get clean

RUN printf "install.packages(c('reshape','ggplot2','wmtsa','parallel'), dep=TRUE, repos='http://cran.us.r-project.org') \n source('http://bioconductor.org/biocLite.R') \n biocLite('rhdf5')" > installpackages.R && \ 
Rscript installpackages.R

RUN git clone https://github.com/lulab/Ribowave.git &&\
	chmod -R +x Ribowave/scripts
