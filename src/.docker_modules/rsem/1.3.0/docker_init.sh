#!/bin/sh
docker pull lbmc/rsem:1.3.0
docker build src/.docker_modules/rsem/1.3.0 -t 'lbmc/rsem:1.3.0'
docker push lbmc/rsem:1.3.0
docker buildx build --platform linux/amd64,linux/arm64 -t "lbmc/rsem:1.3.0" --push src/.docker_modules/rsem/1.3.0
