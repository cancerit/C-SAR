FROM rocker/r-ubuntu:20.04 as builder

# Locale
ENV LC_ALL C
ENV LC_ALL C.UTF-8
ENV LANG C.UTF-8

# Prevent interactive options
ENV DEBIAN_FRONTEND=noninteractive

ENV VER_MAGECK_MAJOR="0.5"
ENV VER_MAGECK_MINOR="9.3"
ENV VER_BAGEL="f9eedca"
ENV VER_NEXTFLOW="21.04.3"

# Run initial system updates
RUN apt-get update && \
  apt-get install -yq --no-install-recommends lsb-release && \
  apt-get update && \
  apt-get install -qy --no-install-recommends \
  apt-transport-https \
  locales \
  libcairo2-dev \
  libcurl4-openssl-dev \
  libxml2-dev \
  libblas-dev \
  libssl-dev \
  libssh2-1-dev \
  r-base \
  r-base-dev \
  python3  \
  python3-dev \
  python3-setuptools \
  python3-pip \
  python3-wheel \
  python3-venv \
  default-jre \
  git \
  curl

ENV OPT /opt/wsi-t113

ENV VIRTUAL_ENV=$OPT/python3
RUN python3 -m venv $VIRTUAL_ENV
ENV PATH="$VIRTUAL_ENV/bin:$OPT/bin:$PATH"

ENV LD_LIBRARY_PATH $OPT/lib
# relying on virtualenv more than this
ENV PYTHONPATH="/usr/src/app:${PYTHONPATH}"
ENV R_LIBS $OPT/R-lib
ENV R_LIBS_USER $R_LIBS

RUN pip3 install --upgrade pip

RUN mkdir -p $R_LIBS_USER

COPY build/install_R_packages.sh build/install_R_packages.sh
RUN bash build/install_R_packages.sh

COPY build/install_R_submodules.sh build/install_R_submodules.sh
COPY submodules/RCRISPR submodules/RCRISPR
RUN bash build/install_R_submodules.sh

COPY build/requirements.txt build/requirements.txt
COPY build/opt-build.sh build/opt-build.sh
RUN bash build/opt-build.sh $OPT

FROM rocker/r-ubuntu:20.04

LABEL maintainer="Victoria Offord <vo1@sanger.ac.uk>" \
      version="1.2.0" \
      description="Nextflow Single CRISPR pipeline container"

RUN apt-get -yq update
RUN apt-get install -yq --no-install-recommends \
  python3 \
  python3-distutils \
  libblas-dev \
  libcurl4 \
  libxml2 \
  libcairo2 \
  r-base \
  default-jre \
  unattended-upgrades && \
  unattended-upgrade -d -v && \
  apt-get remove -yq unattended-upgrades && \
  apt-get autoremove -yq

ENV OPT /opt/wsi-t113

ENV VIRTUAL_ENV=$OPT/python3
ENV PATH="$VIRTUAL_ENV/bin:$OPT/bin:$PATH"

ENV LD_LIBRARY_PATH $OPT/lib
# relying on virtualenv more than this
ENV PYTHONPATH="/usr/src/app:${PYTHONPATH}"
ENV R_LIBS $OPT/R-lib
ENV R_LIBS_USER $R_LIBS

ENV LC_ALL C
ENV LC_ALL C.UTF-8
ENV LANG C.UTF-8
ENV DISPLAY=:0

RUN mkdir -p $OPT
COPY --from=builder $OPT $OPT
RUN mkdir $OPT/c-sar
COPY . $OPT/c-sar

# Copy RCRISPR scripts into pipeline bin directory
RUN ln -s $OPT/c-sar/submodules/RCRISPR/exec/*.R $OPT/bin/

# sort out permissions
RUN find $OPT/c-sar -exec chmod +rx {} \;

RUN adduser --disabled-password --gecos '' ubuntu && chsh -s /bin/bash && mkdir -p /home/ubuntu
USER ubuntu

WORKDIR /home/ubuntu

# ENV NXF_HOME
# ENTRYPOINT ["nextflow", "run", "/opt/wsi-t113/c-sar/main.nf"]

CMD ["/bin/bash"]
