#!/bin/sh
docker pull lbmc/bwa:0.7.17
docker build src/.docker_modules/bwa/0.7.17 -t 'lbmc/bwa:0.7.17'
docker push lbmc/bwa:0.7.17
docker buildx build --platform linux/amd64,linux/arm64 -t "lbmc/bwa:0.7.17" --push src/.docker_modules/bwa/0.7.17
