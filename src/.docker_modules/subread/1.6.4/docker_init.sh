#!/bin/sh
docker pull lbmc/subread:1.6.4
docker build src/.docker_modules/subread/1.6.4 -t 'lbmc/subread:1.6.4'
docker push lbmc/subread:1.6.4
docker buildx build --platform linux/amd64,linux/arm64 -t "lbmc/subread:1.6.4" --push src/.docker_modules/subread/1.6.4
