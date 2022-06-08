#!/bin/sh
docker pull lbmc/music:6613c53
# docker build src/.docker_modules/music/6613c53 -t 'lbmc/music:6613c53'
# docker push lbmc/music:6613c53
docker buildx build --platform linux/amd64,linux/arm64 -t "lbmc/music:6613c53" --push src/.docker_modules/music/6613c53
