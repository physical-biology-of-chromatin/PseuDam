#!/bin/sh
docker pull lbmc/alntools:dd96682
docker build src/.docker_modules/alntools/dd96682 -t 'lbmc/alntools:dd96682'
docker push lbmc/alntools:dd96682
