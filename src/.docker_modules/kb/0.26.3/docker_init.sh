#!/bin/sh
docker pull lbmc/kb:0.26.3
docker build src/.docker_modules/kb/0.26.3 -t 'lbmc/kb:0.26.3'
docker push lbmc/kb:0.26.3
docker buildx build --platform linux/amd64,linux/arm64 -t "lbmc/kb:0.26.3" --push src/.docker_modules/kb/0.26.3
