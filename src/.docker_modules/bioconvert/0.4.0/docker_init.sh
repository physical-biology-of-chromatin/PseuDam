#!/bin/sh
docker pull lbmc/bioconvert:0.4.0
docker build src/.docker_modules/bioconvert/0.4.0 -t 'lbmc/bioconvert:0.4.0'
docker push lbmc/bioconvert:0.4.0
