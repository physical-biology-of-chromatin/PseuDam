#!/bin/sh
docker pull lbmc/rasusa:0.6.0
docker build src/.docker_modules/rasusa/0.6.0 -t 'lbmc/rasusa:0.6.0'
docker push lbmc/rasusa:0.6.0
