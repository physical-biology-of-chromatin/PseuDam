FROM samtools:1.7
MAINTAINER Laurent Modolo

ENV MUSIC_VERSION=6613c53
ENV PACKAGES g++

RUN apk update && \
    apk add ${PACKAGES}

RUN git clone https://github.com/gersteinlab/MUSIC.git && \
  cd MUSIC && \
  git checkout ${MUSIC_VERSION} && \
  make clean && \
  make  && \
  cd .. && \
  mv MUSIC/bin/MUSIC /usr/bin/ && \
  mv MUSIC/bin/generate_multimappability_signal.csh /usr/bin/ && \
  mv MUSIC/bin/run_MUSIC.csh /usr/bin/ && \
  rm -Rf MUSIC

RUN chmod +x /usr/bin/*

