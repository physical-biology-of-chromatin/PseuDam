#!/bin/sh
docker build src/.docker_modules/ucsc/375/ -t 'lbmc/ucsc:375'
docker push lbmc/ucsc:375
