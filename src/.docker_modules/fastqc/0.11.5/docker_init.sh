#!/bin/sh
docker pull lbmc/fastqc:0.11.5
docker build src/.docker_modules/fastqc/0.11.5 -t 'lbmc/fastqc:0.11.5'
docker push lbmc/fastqc:0.11.5
docker buildx build --platform linux/amd64,linux/arm64 -t "lbmc/fastqc:0.11.5" --push src/.docker_modules/fastqc/0.11.5
