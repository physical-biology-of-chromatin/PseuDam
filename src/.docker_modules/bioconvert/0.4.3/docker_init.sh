#!/bin/sh
docker pull lbmc/bioconvert:0.4.3
docker build src/.docker_modules/bioconvert/0.4.3 -t 'lbmc/bioconvert:0.4.3'
docker push lbmc/bioconvert:0.4.3
docker buildx build --platform linux/amd64,linux/arm64 -t "lbmc/bioconvert:0.4.3" --push src/.docker_modules/bioconvert/0.4.3
