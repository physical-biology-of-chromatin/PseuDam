#!/bin/sh
docker pull lbmc/gffread:0.12.2
docker build src/.docker_modules/gffread/0.12.2 -t 'lbmc/gffread:0.12.2'
docker push lbmc/gffread:0.12.2
