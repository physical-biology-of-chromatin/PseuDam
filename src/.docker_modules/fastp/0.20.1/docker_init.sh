#!/bin/sh
docker pull lbmc/fastp:0.20.1
docker build src/.docker_modules/fastp/0.20.1 -t 'lbmc/fastp:0.20.1'
docker push lbmc/fastp:0.20.1
