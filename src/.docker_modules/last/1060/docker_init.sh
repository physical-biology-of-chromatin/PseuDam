#!/bin/sh
docker pull lbmc/last:1060
docker build src/.docker_modules/last/1060/ -t 'lbmc/last:1060'
docker push lbmc/last:1060
docker buildx build --platform linux/amd64,linux/arm64 -t "lbmc/last:1060" --push src/.docker_modules/last/1060
