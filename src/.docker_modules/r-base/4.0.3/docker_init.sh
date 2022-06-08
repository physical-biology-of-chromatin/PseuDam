#!/bin/sh
docker pull lbmc/r-base:4.0.3
# docker build src/.docker_modules/r-base/4.0.3 -t 'lbmc/r-base:4.0.3'
# docker push lbmc/r-base:4.0.3
docker buildx build --platform linux/amd64,linux/arm64 -t "lbmc/r-base:4.0.3" --push src/.docker_modules/r-base/4.0.3
