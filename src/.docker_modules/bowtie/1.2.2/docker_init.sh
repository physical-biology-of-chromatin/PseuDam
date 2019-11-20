#!/bin/sh
docker build src/.docker_modules/bowtie/1.2.2 -t 'lbmc/bowtie:1.2.2'
docker push lbmc/bowtie:1.2.2
