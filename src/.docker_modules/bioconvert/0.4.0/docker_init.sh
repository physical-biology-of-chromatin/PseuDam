#!/bin/sh
docker pull lbmc/bioconvert:0.4.0
docker build src/.docker_modules/bioconvert/0.4.0 -t 'lbmc/bioconvert:0.4.0'
docker push lbmc/bioconvert:0.4.0
docker buildx build --platform linux/amd64,linux/arm64 -t "lbmc/bioconvert:0.4.0" --push src/.docker_modules/bioconvert/0.4.0
