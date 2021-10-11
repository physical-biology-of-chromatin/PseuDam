#!/bin/sh
docker pull lbmc/ucsc:407
docker build src/.docker_modules/ucsc/407/ -t 'lbmc/ucsc:407'
docker push lbmc/ucsc:407
docker buildx build --platform linux/amd64,linux/arm64 -t "lbmc/ucsc:407" --push src/.docker_modules/ucsc/407
