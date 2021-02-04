#!/bin/sh
docker pull lbmc/ucsc:407
docker build src/.docker_modules/ucsc/407/ -t 'lbmc/ucsc:407'
docker push lbmc/ucsc:407
