#!/bin/sh
docker pull lbmc/ribotricer:1.3.2
docker build src/.docker_modules/ribotricer/1.3.2 -t 'lbmc/ribotricer:1.3.2'
docker push lbmc/ribotricer:1.3.2
