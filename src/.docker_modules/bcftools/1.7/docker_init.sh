#!/bin/sh
docker pull lbmc/bcftools:1.7
# docker build src/.docker_modules/bcftools/1.7 -t 'lbmc/bcftools:1.7'
# docker push lbmc/bcftools:1.7
docker buildx build --platform linux/amd64,linux/arm64 -t "lbmc/bcftools:1.7" --push src/.docker_modules/bcftools/1.7
