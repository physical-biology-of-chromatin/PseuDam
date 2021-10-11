#!/bin/sh
docker pull lbmc/flexi_splitter:1.0.2
docker build src/.docker_modules/flexi_splitter/1.0.2 -t 'lbmc/flexi_splitter:1.0.2'
docker push lbmc/flexi_splitter:1.0.2
docker buildx build --platform linux/amd64,linux/arm64 -t "lbmc/flexi_splitter:1.0.2" --push src/.docker_modules/flexi_splitter/1.0.2
