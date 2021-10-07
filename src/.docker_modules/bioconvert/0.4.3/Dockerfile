FROM python:3.9-slim 

ENV BIOCONVERT_VERSION="0.4.3"

RUN apt update && apt install -y procps gcc zlib1g-dev && pip3 install bioconvert==${BIOCONVERT_VERSION}

CMD [ "bash" ]
