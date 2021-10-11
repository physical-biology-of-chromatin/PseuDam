#!/bin/sh
docker pull lbmc/ucsc:375
docker build src/.docker_modules/ucsc/375/ -t 'lbmc/ucsc:375'
docker push lbmc/ucsc:375
docker buildx build --platform linux/amd64,linux/arm64 -t "lbmc/ucsc:375" --push src/.docker_modules/ucsc/375
