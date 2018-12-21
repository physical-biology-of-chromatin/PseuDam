FROM alpine:3.8
MAINTAINER Laurent Modolo

ENV UCSC_VERSION=375
ENV PACKAGES rsync \
             bash

RUN apk update && \
    apk add ${PACKAGES}

RUN rsync -aP rsync://hgdownload.soe.ucsc.edu/genome/admin/exe/linux.x86_64/ /usr/bin/

