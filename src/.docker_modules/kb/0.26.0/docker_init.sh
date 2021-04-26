#!/bin/sh
docker pull lbmc/kb:0.26.0
docker build src/.docker_modules/kb/0.26.0 -t 'lbmc/kb:0.26.0'
docker push lbmc/kb:0.26.0
