#!/bin/sh
docker pull lbmc/pandoc:2.11
# docker build src/.docker_modules/pandoc/2.11 -t 'lbmc/pandoc:2.11'
# docker push lbmc/pandoc:2.11
docker buildx build --platform linux/amd64,linux/arm64 -t "lbmc/pandoc:2.11" --push src/.docker_modules/pandoc/2.11
