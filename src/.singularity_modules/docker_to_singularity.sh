# install spython
# sudo pip3 install spython

find src/docker_modules/ -name "Dockerfile" | \
  perl -pe "s/docker/singularity/g" | \
  perl -pe "s/Dockerfile//g" | \
  awk '{system("mkdir -p " $0)}'

find src/docker_modules/ -name "Dockerfile" | \
  perl -pe "s/(^.*$)/spython recipe \1 > \1/g" | \
  perl -pe "s/(^.*)docker_modules(.*$)/\1singularity_modules\2/g" | \
  perl -pe "s/(^.*\/([^\/]*)\/[^\/]*\/)Dockerfile$/\1\2.def/g" | \
  awk '{system($0)}'
