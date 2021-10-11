#!/bin/sh
docker pull lbmc/kallisto:0.43.1
docker build src/.docker_modules/kallisto/0.43.1 -t 'lbmc/kallisto:0.43.1'
docker push lbmc/kallisto:0.43.1
docker buildx build --platform linux/amd64,linux/arm64 -t "lbmc/kallisto:0.43.1" --push src/.docker_modules/kallisto/0.43.1
