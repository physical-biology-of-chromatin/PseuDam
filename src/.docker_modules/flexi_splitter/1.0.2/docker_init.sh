#!/bin/sh
docker pull lbmc/flexi_splitter:1.0.2
docker build src/.docker_modules/flexi_splitter/1.0.2 -t 'lbmc/flexi_splitter:1.0.2'
docker push lbmc/flexi_splitter:1.0.2
