#!/bin/sh
docker pull lbmc/r-base:4.0.0
docker build src/.docker_modules/r-base/4.0.0 -t 'lbmc/r-base:4.0.0'
docker push lbmc/r-base:4.0.0
docker buildx build --platform linux/amd64,linux/arm64 -t "lbmc/r-base:4.0.0" --push src/.docker_modules/r-base/4.0.0
