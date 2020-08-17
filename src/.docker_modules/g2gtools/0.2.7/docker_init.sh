#!/bin/sh
docker pull lbmc/g2gtools:0.2.7
docker build src/.docker_modules/g2gtools/0.2.7 -t 'lbmc/g2gtools:0.2.7'
docker push lbmc/g2gtools:0.2.7
