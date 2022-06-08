#!/bin/sh
docker pull lbmc/canu:1.6
# docker build src/.docker_modules/canu/1.6 -t 'lbmc/canu:1.6'
# docker push lbmc/canu:1.6
docker buildx build --platform linux/amd64,linux/arm64 -t "lbmc/canu:1.6" --push src/.docker_modules/canu/1.6
