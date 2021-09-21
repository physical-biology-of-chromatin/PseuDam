#!/bin/sh
docker pull lbmc/danpos3:2f7f223
docker build src/.docker_modules/danpos3/2f7f223 -t 'lbmc/danpos3:2f7f223'
docker push lbmc/danpos3:2f7f223
