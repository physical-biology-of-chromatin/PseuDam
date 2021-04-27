#!/bin/sh
docker pull lbmc/deeptools:3.5.1
docker build src/.docker_modules/deeptools/3.5.1 -t 'lbmc/deeptools:3.5.1'
docker push lbmc/deeptools:3.5.1
