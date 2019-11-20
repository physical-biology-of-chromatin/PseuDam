#!/bin/sh
docker build src/.docker_modules/sambamba/0.6.9 -t 'lbmc/sambamba:0.6.9'
docker push lbmc/sambamba:0.6.9
docker push lbmc/sambamba:0.6.9
