#!/bin/sh
docker pull lbmc/crossmap:0.4.1
# docker build src/.docker_modules/crossmap/0.4.1/ -t 'lbmc/crossmap:0.4.1'
# docker push lbmc/crossmap:0.4.1
docker buildx build --platform linux/amd64,linux/arm64 -t "lbmc/crossmap:0.4.1" --push src/.docker_modules/crossmap/0.4.1
