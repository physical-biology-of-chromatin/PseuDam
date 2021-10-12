#!/bin/sh
docker pull lbmc/deeptools:3.0.2
# docker build src/.docker_modules/deeptools/3.0.2 -t 'lbmc/deeptools:3.0.2'
# docker push lbmc/deeptools:3.0.2
docker buildx build --platform linux/amd64,linux/arm64 -t "lbmc/deeptools:3.0.2" --push src/.docker_modules/deeptools/3.0.2
