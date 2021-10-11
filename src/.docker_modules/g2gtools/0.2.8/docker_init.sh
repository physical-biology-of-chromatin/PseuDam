#!/bin/sh
docker pull lbmc/g2gtools:0.2.8
docker build src/.docker_modules/g2gtools/0.2.8 -t 'lbmc/g2gtools:0.2.8'
docker push lbmc/g2gtools:0.2.8
docker buildx build --platform linux/amd64,linux/arm64 -t "lbmc/g2gtools:0.2.8" --push src/.docker_modules/g2gtools/0.2.8
