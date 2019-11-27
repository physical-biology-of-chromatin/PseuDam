#!/bin/sh
docker pull lbmc/htseq:0.11.2
docker build src/.docker_modules/htseq/0.11.2 -t 'lbmc/htseq:0.11.2'
docker push lbmc/htseq:0.11.2
