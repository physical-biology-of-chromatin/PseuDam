#!/bin/sh
docker pull lbmc/bowtie2:2.3.4.1
docker build src/.docker_modules/bowtie2/2.3.4.1 -t 'lbmc/bowtie2:2.3.4.1'
docker push lbmc/bowtie2:2.3.4.1
docker buildx build --platform linux/amd64,linux/arm64 -t "lbmc/bowtie2:2.3.4.1" --push src/.docker_modules/bowtie2/2.3.4.1
