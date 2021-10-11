#!/bin/sh
docker pull lbmc/emase-zero:0.3.1
docker build src/.docker_modules/emase-zero/0.3.1 -t 'lbmc/emase-zero:0.3.1'
docker push lbmc/emase-zero:0.3.1
docker buildx build --platform linux/amd64,linux/arm64 -t "lbmc/emase-zero:0.3.1" --push src/.docker_modules/emase-zero/0.3.1
