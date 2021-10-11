#!/bin/sh
docker pull lbmc/deeptools:3.5.0
docker build src/.docker_modules/deeptools/3.5.0 -t 'lbmc/deeptools:3.5.0'
docker push lbmc/deeptools:3.5.0
docker buildx build --platform linux/amd64,linux/arm64 -t "lbmc/deeptools:3.5.0" --push src/.docker_modules/deeptools/3.5.0
