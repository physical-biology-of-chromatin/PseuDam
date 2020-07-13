#!/bin/sh
docker pull lbmc/crossmap:0.4.1
docker build src/.docker_modules/crossmap/0.4.1/ -t 'lbmc/crossmap:0.4.1'
docker push lbmc/crossmap:0.4.1
