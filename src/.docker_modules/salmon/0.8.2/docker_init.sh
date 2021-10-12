#!/bin/sh
docker pull lbmc/salmon:0.8.2
# docker build src/.docker_modules/salmon/0.8.2 -t 'lbmc/salmon:0.8.2'
# docker push lbmc/salmon:0.8.2
docker buildx build --platform linux/amd64,linux/arm64 -t "lbmc/salmon:0.8.2" --push src/.docker_modules/salmon/0.8.2
