FROM debian:stretch
MAINTAINER Rémi SERAPHIN

ENV UMI_TOOLS_VERSION=1.0.0
ENV PACKAGES="bash \
             python3 \
             python3-dev \
             python3-pip \
             procps \
             gcc"

RUN apt-get update && \
apt-get install -y --no-install-recommends ${PACKAGES} && \
apt-get clean

RUN pip3 install setuptools
RUN pip3 install umi_tools==${UMI_TOOLS_VERSION}

