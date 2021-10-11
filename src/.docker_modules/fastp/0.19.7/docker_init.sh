#!/bin/sh
docker pull lbmc/fastp:0.19.7
docker build src/.docker_modules/fastp/0.19.7 -t 'lbmc/fastp:0.19.7'
docker push lbmc/fastp:0.19.7
docker buildx build --platform linux/amd64,linux/arm64 -t "lbmc/fastp:0.19.7" --push src/.docker_modules/fastp/0.19.7
