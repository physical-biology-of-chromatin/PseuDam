#!/bin/sh
docker pull lbmc/r-base:3.6.2
docker build src/.docker_modules/r-base/3.6.2 -t 'lbmc/r-base:3.6.2'
docker push lbmc/r-base:3.6.2
docker buildx build --platform linux/amd64,linux/arm64 -t "lbmc/r-base:3.6.2" --push src/.docker_modules/r-base/3.6.2
