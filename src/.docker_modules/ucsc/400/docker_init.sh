#!/bin/sh
docker pull lbmc/ucsc:400
# docker build src/.docker_modules/ucsc/400/ -t 'lbmc/ucsc:400'
# docker push lbmc/ucsc:400
docker buildx build --platform linux/amd64,linux/arm64 -t "lbmc/ucsc:400" --push src/.docker_modules/ucsc/400
