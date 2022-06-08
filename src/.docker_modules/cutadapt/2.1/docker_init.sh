#!/bin/sh
docker pull lbmc/cutadapt:2.1
# docker build src/.docker_modules/cutadapt/2.1 -t 'lbmc/cutadapt:2.1'
# docker push lbmc/cutadapt:2.1
docker buildx build --platform linux/amd64,linux/arm64 -t "lbmc/cutadapt:2.1" --push src/.docker_modules/cutadapt/2.1
