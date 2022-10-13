FROM rocker/tidyverse:4.1.1

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update \
 && apt-get -y install \
    build-essential \
    libcairo2-dev \
    libcurl4-gnutls-dev \
    libncurses-dev \
    libgeos-dev \
    libglpk-dev \
    libgsl-dev \
    libmagick++-dev \
    libxml2-dev \
    libssl-dev \
    python3-dev \
    python3-pip \
 && rm -rf /var/lib/apt/lists

COPY install_R_packages.R /opt

RUN Rscript /opt/install_R_packages.R \
 && rm /opt/install_R_packages.R

COPY bin /opt

CMD ["/bin/bash"]