#!/bin/sh
docker pull lbmc/bioawk:1.0
docker build src/.docker_modules/bioawk/1.0 -t 'lbmc/bioawk:1.0'
docker push lbmc/bioawk:1.0
docker buildx build --platform linux/amd64,linux/arm64 -t "lbmc/bioawk:1.0" --push src/.docker_modules/bioawk/1.0
