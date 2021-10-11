#!/bin/sh
docker pull lbmc/gatk:3.8.0
docker build src/.docker_modules/gatk/3.8.0 -t 'lbmc/gatk:3.8.0'
docker push lbmc/gatk:3.8.0
docker buildx build --platform linux/amd64,linux/arm64 -t "lbmc/gatk:3.8.0" --push src/.docker_modules/gatk/3.8.0
