FROM python:3.8-alpine

RUN apk update \
    && apk add make automake gcc g++ subversion python3-dev

RUN pip install numpy\
    && pip install biopython

COPY script/ /script

RUN adduser -D -u 1000 finder
RUN chown -R finder /script
USER finder

ENTRYPOINT [ "python","script/GATC_finder.py" ]