FROM ubuntu:18.04
MAINTAINER Laurent Modolo

ENV BOWTIE2_VERSION=2.3.4.1
ENV PACKAGES bowtie2=${BOWTIE2_VERSION}*

RUN apt-get update && \
    apt-get install -y --no-install-recommends ${PACKAGES} && \
    apt-get clean
