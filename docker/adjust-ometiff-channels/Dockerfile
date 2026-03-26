FROM ubuntu:24.04

ENV DEBIAN_FRONTEND noninteractive

RUN apt-get update \
 && apt-get -y install \
    git \
    python3 \
    python3-pip \
    python3-venv \
 && rm -rf /var/cache/apt/*

COPY docker/adjust-ometiff-channels/requirements.txt /opt
RUN python3 -m venv /opt/venv \
 && . /opt/venv/bin/activate \
 && pip install -r /opt/requirements.txt \
 && rm -rf /root/.cache/pip /opt/requirements.txt

COPY bin/common /opt
COPY bin/convert_to_ometiff /opt

env PATH /opt/venv/bin:${PATH}
