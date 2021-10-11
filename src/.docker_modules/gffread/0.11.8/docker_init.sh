#!/bin/sh
docker pull lbmc/gffread:0.11.8
docker build src/.docker_modules/gffread/0.11.8 -t 'lbmc/gffread:0.11.8'
docker push lbmc/gffread:0.11.8
docker buildx build --platform linux/amd64,linux/arm64 -t "lbmc/gffread:0.11.8" --push src/.docker_modules/gffread/0.11.8
