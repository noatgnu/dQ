FROM python:3.9

EXPOSE 8000

RUN mkdir /app
RUN mkdir /app/src
RUN mkdir /app/data

WORKDIR /app/src

RUN apk update
RUN apk add git

RUN git clone https://github.com/noatgnu/dQ.git

WORKDIR /app/src/dQ

RUN pip -r install requirements.txt

ENTRYPOINT ["python", "main.py"]
