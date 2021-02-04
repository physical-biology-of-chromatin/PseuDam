#!/bin/sh
docker pull lbmc/deeptools:3.5.0
docker build src/.docker_modules/deeptools/3.5.0 -t 'lbmc/deeptools:3.5.0'
docker push lbmc/deeptools:3.5.0
