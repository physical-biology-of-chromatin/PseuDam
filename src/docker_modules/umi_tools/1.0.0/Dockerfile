FROM debian:stretch
MAINTAINER Rémi SERAPHIN

ENV UMI_TOOLS_VERSION=1.0.0
ENV PACKAGES bash \
             python3 \
             python3-pip \
             gcc

RUN apt update && \
    apt install ${PACKAGES} -y

RUN pip3 install umi_tools==${UMI_TOOLS_VERSION}
