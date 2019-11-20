#!/bin/sh
docker build src/.docker_modules/subread/1.6.4 -t 'lbmc/subread:1.6.4'
docker push lbmc/subread:1.6.4
