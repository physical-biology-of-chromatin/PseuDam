#!/bin/sh
docker pull lbmc/freebayes:1.3.2
docker build src/.docker_modules/freebayes/1.3.2/ -t 'lbmc/freebayes:1.3.2'
docker push lbmc/freebayes:1.3.2
docker buildx build --platform linux/amd64,linux/arm64 -t "lbmc/freebayes:1.3.2" --push src/.docker_modules/freebayes/1.3.2
