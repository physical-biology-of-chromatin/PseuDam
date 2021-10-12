#!/bin/sh
docker pull lbmc/kallisto:0.44.0
# docker build src/.docker_modules/kallisto/0.44.0 -t 'lbmc/kallisto:0.44.0'
# docker push lbmc/kallisto:0.44.0
docker buildx build --platform linux/amd64,linux/arm64 -t "lbmc/kallisto:0.44.0" --push src/.docker_modules/kallisto/0.44.0
