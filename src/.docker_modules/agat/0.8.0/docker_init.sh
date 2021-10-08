#!/bin/sh
docker pull lbmc/agat:0.8.0
docker build src/.docker_modules/agat/0.8.0 -t 'lbmc/agat:0.8.0'
docker push lbmc/agat:0.8.0
