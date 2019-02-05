FROM debian:jessie
MAINTAINER Laurent Modolo

ENV PACKAGES apt-utils \
             curl \
             build-essential

RUN apt-get update && \
    apt-get install -y ${PACKAGES}

ENV UCSC_VERSION=375

RUN curl -k -L http://hgdownload.soe.ucsc.edu/admin/exe/userApps.v${UCSC_VERSION}.src.tgz -o userApps.v${UCSC_VERSION}.src.tgz &&\
tar xvf userApps.v${UCSC_VERSION}.src.tgz &&\
cd userApps/ && \
make
