#!/bin/sh
docker pull lbmc/hisat2:2.0.0
# docker build src/.docker_modules/hisat2/2.0.0 -t 'lbmc/hisat2:2.0.0'
# docker push lbmc/hisat2:2.0.0
docker buildx build --platform linux/amd64,linux/arm64 -t "lbmc/hisat2:2.0.0" --push src/.docker_modules/hisat2/2.0.0
