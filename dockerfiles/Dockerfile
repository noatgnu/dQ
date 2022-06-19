FROM python:3.9

ARG DOCKER_IMAGE_VERSION=unknown

LABEL \
      org.label-schema.name="dQ" \
      org.label-schema.description="Docker container for dQ" \
      org.label-schema.version="$DOCKER_IMAGE_VERSION" \
      org.label-schema.vcs-url="https://github.com/noatgnu/dQ" \
      org.label-schema.schema-version="1.0"

EXPOSE 8000

RUN mkdir /app
RUN mkdir /app/src
RUN mkdir /app/data

WORKDIR /app/src

RUN apt-get update
RUN apt-get upgrade -y
RUN apt-get install git

RUN git clone https://github.com/noatgnu/dQ.git

WORKDIR /app/src/dQ

RUN pip install -r requirements.txt

VOLUME ["/app/data"]

ENTRYPOINT ["python", "main.py"]
