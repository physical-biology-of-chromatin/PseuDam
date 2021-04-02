#!/bin/sh
docker pull lbmc/macs3:3.0.0a6
docker build src/.docker_modules/macs3/3.0.0a6 -t 'lbmc/macs3:3.0.0a6'
docker push lbmc/macs3:3.0.0a6
