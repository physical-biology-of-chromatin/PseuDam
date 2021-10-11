#!/bin/sh
docker pull lbmc/samblaster:0.1.24
docker build src/.docker_modules/samblaster/0.1.24 -t 'lbmc/samblaster:0.1.24'
docker push lbmc/samblaster:0.1.24
docker buildx build --platform linux/amd64,linux/arm64 -t "lbmc/samblaster:0.1.24" --push src/.docker_modules/samblaster/0.1.24
