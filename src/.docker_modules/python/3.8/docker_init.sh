#!/bin/sh
docker pull lbmc/python:3.8
docker build src/.docker_modules/python/3.8 -t 'lbmc/python:3.8'
docker push lbmc/python:3.8
docker buildx build --platform linux/amd64,linux/arm64 -t "lbmc/python:3.8" --push src/.docker_modules/python/3.8
