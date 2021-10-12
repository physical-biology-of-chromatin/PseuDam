#!/bin/sh
docker pull lbmc/fastp:0.20.1
# docker build src/.docker_modules/fastp/0.20.1 -t 'lbmc/fastp:0.20.1'
# docker push lbmc/fastp:0.20.1
docker buildx build --platform linux/amd64,linux/arm64 -t "lbmc/fastp:0.20.1" --push src/.docker_modules/fastp/0.20.1
