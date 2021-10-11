#!/bin/sh
docker pull lbmc/sambamba:0.6.9
docker build src/.docker_modules/sambamba/0.6.9 -t 'lbmc/sambamba:0.6.9'
docker push lbmc/sambamba:0.6.9
docker buildx build --platform linux/amd64,linux/arm64 -t "lbmc/sambamba:0.6.9" --push src/.docker_modules/sambamba/0.6.9
