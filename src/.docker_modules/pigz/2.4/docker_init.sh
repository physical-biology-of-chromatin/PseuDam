#!/bin/sh
docker pull lbmc/pigz:2.4
# docker build src/.docker_modules/pigz/2.4 -t 'lbmc/pigz:2.4'
# docker push lbmc/pigz:2.4
docker buildx build --platform linux/amd64,linux/arm64 -t "lbmc/pigz:2.4" --push src/.docker_modules/pigz/2.4
