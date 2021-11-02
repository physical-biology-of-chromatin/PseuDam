#!/bin/sh
docker pull lbmc/multiqc:1.11
# docker build src/.docker_modules/multiqc/1.11 -t 'lbmc/multiqc:1.11'
# docker push lbmc/multiqc:1.11
docker buildx build --platform linux/amd64,linux/arm64 -t "lbmc/multiqc:1.11" --push src/.docker_modules/multiqc/1.11
