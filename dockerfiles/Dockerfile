FROM python:3.9

ARG DOCKER_IMAGE_VERSION=0.1

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
RUN apt-get install git redis -y
RUN apt-get install -y --no-install-recommends r-base r-base-dev && rm -rf /var/lib/apt/lists/*
RUN R -e "install.packages('BiocManager',dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "library('BiocManager');BiocManager::install('limma')"

RUN git clone https://github.com/noatgnu/dQ.git

WORKDIR /app/src/dQ

RUN pip install --no-cache-dir -r requirements.txt

VOLUME ["/app/data"]

