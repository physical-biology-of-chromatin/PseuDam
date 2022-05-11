FROM python:3.8-alpine

RUN apk update \
    && apk add g++ \
    && apk add --no-cache bash

RUN pip install numpy\
    && pip install biopython

COPY gatc_finder.py /usr/bin/

RUN chmod +x /usr/bin/*.py


