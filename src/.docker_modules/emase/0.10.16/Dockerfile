FROM python:2.7
MAINTAINER Laurent Modolo

ENV EMASE_VERSION=0.10.16

RUN apt-get update \
    && apt-get install -y python-h5py \
    && pip install numpy==1.8.2 \
    && pip install numexpr==2.3.1 \
    && pip install cython>=0.13 \
    && pip install tables==3.1.0 \
    && pip install emase==${EMASE_VERSION}

CMD ["bash"]

