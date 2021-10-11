#!/bin/sh
docker pull lbmc/rasusa:0.6.0
docker build src/.docker_modules/rasusa/0.6.0 -t 'lbmc/rasusa:0.6.0'
docker push lbmc/rasusa:0.6.0
docker buildx build --platform linux/amd64,linux/arm64 -t "lbmc/rasusa:0.6.0" --push src/.docker_modules/rasusa/0.6.0
docker buildx build --platform linux/amd64,linux/arm64 -t 'lbmc/rasusa:0.6.0' --push src/.docker_modules/rasusa/0.6.0
