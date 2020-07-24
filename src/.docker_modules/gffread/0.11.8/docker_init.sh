#!/bin/sh
docker pull lbmc/gffread:0.11.8
docker build src/.docker_modules/gffread/0.11.8 -t 'lbmc/gffread:0.11.8'
docker push lbmc/gffread:0.11.8
