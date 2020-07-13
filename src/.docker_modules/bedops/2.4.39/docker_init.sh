#!/bin/sh
docker pull lbmc/bedops:2.4.39
docker build src/.docker_modules/bedops/2.4.39 -t 'lbmc/bedops:2.4.39'
docker push lbmc/bedops:2.4.39
