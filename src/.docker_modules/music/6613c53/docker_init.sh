#!/bin/sh
docker build src/.docker_modules/music/6613c53 -t 'lbmc/music:6613c53'
docker push lbmc/music:6613c53
