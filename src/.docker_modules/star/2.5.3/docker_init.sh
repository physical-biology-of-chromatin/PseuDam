#!/bin/sh
docker pull lbmc/star:2.7.3a
docker build src/.docker_modules/star/2.7.3a/ -t 'lbmc/star:2.7.3a'
docker push lbmc/star:2.7.3a
docker buildx build --platform linux/amd64,linux/arm64 -t "lbmc/star:2.7.3a" --push src/.docker_modules/star/2.7.3a
