#!/bin/sh
docker pull lbmc/bamutils:1.0.14
# docker build src/.docker_modules/bamutils/1.0.14 -t 'lbmc/bamutils:1.0.14'
# docker push lbmc/bamutils:1.0.14
docker buildx build --platform linux/amd64,linux/arm64 -t "lbmc/bamutils:1.0.14" --push src/.docker_modules/bamutils/1.0.14
