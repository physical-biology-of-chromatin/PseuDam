#!/bin/sh
docker pull lbmc/emase:0.10.16
docker build src/.docker_modules/emase/0.10.16 -t 'lbmc/emase:0.10.16'
docker push lbmc/emase:0.10.16
