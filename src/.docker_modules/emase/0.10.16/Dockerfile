FROM python:2.7
MAINTAINER Laurent Modolo

ENV EMASE_VERSION=0.10.16

RUN 
apt-get update \
    && apt-get install -y python-h5py libhdf5-dev \
    && pip install tables \
    && pip install scipy \
    && pip install pysam \
    && pip install numpy==1.8.2 \
    && pip install numexpr==2.3.1 \
    && pip install cython>=0.13 \
    && pip install tables==3.1.0 \
    && pip install emase==${EMASE_VERSION}

pip install numpy==1.8.2
pip install numexpr==2.3.1 
pip install scipy==0.13.3 
pip install pysam>=0.6 
pip install cython>=0.13 
pip install tables
pip install biopython>=1.63

CMD ["bash"]

