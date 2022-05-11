FROM perl:5.34

MAINTAINER Nathan Lecouvreur

RUN apk update \
    && apk add --no-cache bash


COPY fastq_splitter.pl /usr/bin/

RUN chmod +x /usr/bin/*.pl