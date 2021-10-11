#!/bin/sh
docker pull lbmc/kallistobustools:0.39.3
docker build src/.docker_modules/kallistobustools/0.39.3 -t 'lbmc/kallistobustools:0.39.3'
docker push lbmc/kallistobustools:0.39.3
docker buildx build --platform linux/amd64,linux/arm64 -t "lbmc/kallistobustools:0.39.3" --push src/.docker_modules/kallistobustools/0.39.3
