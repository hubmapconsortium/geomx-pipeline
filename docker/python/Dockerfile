FROM ubuntu:24.04
ENV DEBIAN_FRONTEND=noninteractive

## System packages
# I'm used to using `apt` but it complains that it doesn't have a stable
#   command-line interface, so let's use `apt-get` instead
RUN apt-get update \
 && apt-get -y install \
    git \
    python3 \
    python3-dev \
    python3-pip \
    python3-venv \
 && rm -rf /var/cache/apt/*

## Python packages
COPY docker/python/requirements.txt /opt/requirements.txt
RUN python3 -m venv /opt/venv \
 && . /opt/venv/bin/activate \
 && pip install -r /opt/requirements.txt \
 && rm -rf /root/.cache/pip /opt/requirements.txt

RUN mkdir /root/.ssh

COPY bin /opt
COPY data /opt

WORKDIR /opt

env PATH /opt/venv/bin:${PATH}
