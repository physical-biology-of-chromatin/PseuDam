#!/bin/sh
docker pull lbmc/umi_tools:1.0.0
docker build src/.docker_modules/umi_tools/1.0.0/ -t 'lbmc/umi_tools:1.0.0'
docker push lbmc/umi_tools:1.0.0
docker buildx build --platform linux/amd64,linux/arm64 -t "lbmc/umi_tools:1.0.0" --push src/.docker_modules/umi_tools/1.0.0
