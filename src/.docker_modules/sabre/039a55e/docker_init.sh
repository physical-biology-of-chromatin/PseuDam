#!/bin/sh
docker pull lbmc/sabre:039a55e
docker build src/.docker_modules/sabre/039a55e -t 'lbmc/sabre:039a55e'
docker push lbmc/sabre:039a55e
