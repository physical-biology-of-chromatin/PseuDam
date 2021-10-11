#!/bin/sh
docker pull lbmc/sratoolkit:2.8.2
docker build src/.docker_modules/sratoolkit/2.8.2 -t 'lbmc/sratoolkit:2.8.2'
docker push lbmc/sratoolkit:2.8.2
docker buildx build --platform linux/amd64,linux/arm64 -t "lbmc/sratoolkit:2.8.2" --push src/.docker_modules/sratoolkit/2.8.2
