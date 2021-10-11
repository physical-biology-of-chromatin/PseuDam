#!/bin/sh
docker pull lbmc/liftover:357
docker build src/.docker_modules/liftover/357/ -t 'lbmc/liftover:357'
docker push lbmc/liftover:357
docker buildx build --platform linux/amd64,linux/arm64 -t "lbmc/liftover:357" --push src/.docker_modules/liftover/357
