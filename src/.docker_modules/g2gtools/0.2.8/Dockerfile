FROM python:3.8-alpine
MAINTAINER Laurent Modolo

ENV G2GTOOLS_VERSION=0.2.8

RUN apk add --update --no-cache bash musl-dev linux-headers g++ cmake make build-base bzip2-dev zlib-dev xz-dev autoconf \
 && wget https://github.com/churchill-lab/g2gtools/archive/v${G2GTOOLS_VERSION}.tar.gz \
 && tar -xvf v${G2GTOOLS_VERSION}.tar.gz \
 && cd g2gtools-${G2GTOOLS_VERSION} \
 && pip install numpy \
 && pip install pysam \
 && make install 


CMD ["bash"]
