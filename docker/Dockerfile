FROM Ubuntu:20.04
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
 && rm -rf /var/cache/apt/*

## Python packages
COPY docker/requirements.txt /opt/requirements.txt
RUN pip3 install --upgrade pip

RUN pip3 install -r /opt/requirements.txt \
 && rm -rf /root/.cache/pip

RUN mkdir /root/.ssh

COPY bin /opt

WORKDIR /opt
