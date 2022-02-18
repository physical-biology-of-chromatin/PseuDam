FROM python:3.8-alpine

RUN apk update \
    && apk add g++

RUN pip install numpy\
    && pip install biopython

COPY gatc_finder.py .

RUN echo > sites.bed



ENTRYPOINT [ "" ]
