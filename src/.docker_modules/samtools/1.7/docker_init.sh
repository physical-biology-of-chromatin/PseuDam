#!/bin/sh
docker pull lbmc/samtools:1.7
docker build src/.docker_modules/samtools/1.7 -t 'lbmc/samtools:1.7'
docker push lbmc/samtools:1.7
docker buildx build --platform linux/amd64,linux/arm64 -t "lbmc/samtools:1.7" --push src/.docker_modules/samtools/1.7
