#!/bin/sh
docker pull lbmc/tophat:2.1.1
# docker build src/.docker_modules/tophat/2.1.1 -t 'lbmc/tophat:2.1.1'
# docker push lbmc/tophat:2.1.1
docker buildx build --platform linux/amd64,linux/arm64 -t "lbmc/tophat:2.1.1" --push src/.docker_modules/tophat/2.1.1
