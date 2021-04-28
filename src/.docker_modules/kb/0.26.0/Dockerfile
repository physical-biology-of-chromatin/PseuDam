FROM python:3.9-slim 

ENV KB_VERSION="0.26.0"

RUN apt update && apt install -y procps && pip3 install kb-python==${KB_VERSION}

CMD [ "bash" ]
