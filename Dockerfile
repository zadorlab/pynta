FROM ubuntu:18.04

LABEL Maciej Gierada 'maciej.gierada@gmail.com'

RUN apt-get update -y && \
    apt-get install -y python3-pip python3-dev

COPY ./requirements.txt /requirements.txt

WORKDIR /

RUN python3 setup.py install --user

COPY . /

ENTRYPOINT [ "python3" ]