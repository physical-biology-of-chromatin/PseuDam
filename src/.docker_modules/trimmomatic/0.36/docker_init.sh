#!/bin/sh
docker build src/.docker_modules/trimmomatic/0.36 -t 'lbmc/trimmomatic:0.36'
docker push lbmc/trimmomatic:0.36
