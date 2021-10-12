#!/bin/sh
docker pull lbmc/bowtie:1.2.2
# docker build src/.docker_modules/bowtie/1.2.2 -t 'lbmc/bowtie:1.2.2'
# docker push lbmc/bowtie:1.2.2
docker buildx build --platform linux/amd64,linux/arm64 -t "lbmc/bowtie:1.2.2" --push src/.docker_modules/bowtie/1.2.2
