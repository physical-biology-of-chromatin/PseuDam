#!/bin/sh
docker pull lbmc/samtools:1.14
# docker build src/.docker_modules/samtools/1.14 -t 'lbmc/samtools:1.14'
# docker push lbmc/samtools:1.14
docker buildx build --platform linux/amd64,linux/arm64 -t "lbmc/samtools:1.14" --push src/.docker_modules/samtools/1.14
