#!/bin/sh
docker pull lbmc/file_handle:0.1.1
docker build src/.docker_modules/file_handle/0.1.1 -t 'lbmc/file_handle:0.1.1'
docker push lbmc/file_handle:0.1.1
docker buildx build --platform linux/amd64,linux/arm64 -t "lbmc/file_handle:0.1.1" --push src/.docker_modules/file_handle/0.1.1
