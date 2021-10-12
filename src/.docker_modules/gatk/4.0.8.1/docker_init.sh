#!/bin/sh
docker pull lbmc/gatk:4.0.8.1
# docker build src/.docker_modules/gatk/4.0.8.1 -t 'lbmc/gatk:4.0.8.1'
# docker push lbmc/gatk:4.0.8.1
docker buildx build --platform linux/amd64,linux/arm64 -t "lbmc/gatk:4.0.8.1" --push src/.docker_modules/gatk/4.0.8.1
