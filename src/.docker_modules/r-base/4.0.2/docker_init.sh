#!/bin/sh
docker pull lbmc/r-base:4.0.2
docker build src/.docker_modules/r-base/4.0.2 -t 'lbmc/r-base:4.0.2'
docker push lbmc/r-base:4.0.2
