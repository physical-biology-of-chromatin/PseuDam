#!/bin/sh
docker pull lbmc/minimap2:2.17
# docker build src/.docker_modules/minimap2/2.17 -t 'lbmc/minimap2:2.17'
# docker push lbmc/minimap2:2.17
docker buildx build --platform linux/amd64,linux/arm64 -t "lbmc/minimap2:2.17" --push src/.docker_modules/minimap2/2.17
