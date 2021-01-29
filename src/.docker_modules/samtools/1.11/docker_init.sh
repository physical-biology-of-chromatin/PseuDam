#!/bin/sh
docker pull lbmc/samtools:1.11
docker build src/.docker_modules/samtools/1.11 -t 'lbmc/samtools:1.11'
docker push lbmc/samtools:1.11
