#!/bin/sh
docker pull lbmc/r-base:4.0.3
docker build src/.docker_modules/r-base/4.0.3 -t 'lbmc/r-base:4.0.3'
docker push lbmc/r-base:4.0.3
