#!/bin/sh
docker pull lbmc/trimmomatic:0.36
# docker build src/.docker_modules/trimmomatic/0.36 -t 'lbmc/trimmomatic:0.36'
# docker push lbmc/trimmomatic:0.36
docker buildx build --platform linux/amd64,linux/arm64 -t "lbmc/trimmomatic:0.36" --push src/.docker_modules/trimmomatic/0.36
