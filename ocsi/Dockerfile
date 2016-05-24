# run as:
#   docker build -t ocsi .

# base everything on a recent Ubuntu
FROM debian:latest

# get system packages up to date then install a basic scientific python
RUN apt-get update && apt-get -y upgrade && \
    apt-get -y install python \
         python-numpy python-scipy python-pandas ttf-bitstream-vera

# add code
RUN mkdir /scripts
COPY scripts /scripts

MAINTAINER Krzysztof Polanski <k.t.polanski@warwick.ac.uk>

# this is where we start
ENTRYPOINT ["bash", "/scripts/ocsi_tarwrapper.sh"]

# if nothing else is specified in the docker call, just run --help
CMD ["--help"]