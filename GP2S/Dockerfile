# run as:
#   docker build -t gp2s .

# base everything on a recent Ubuntu
FROM debian:latest

# get system packages up to date then install a basic scientific python
RUN apt-get update && apt-get -y upgrade && \
    apt-get -y install python python-pip python-dev build-essential \
    python-numpy python-scipy python-matplotlib ttf-bitstream-vera

# upgrade pandas
RUN pip install --upgrade pandas
# downgrade numpy, as the pandas upgraded version somehow breaks time local mode
RUN pip install -I numpy==1.8.2

# add code
RUN mkdir /scripts
COPY scripts /scripts

MAINTAINER Krzysztof Polanski <k.t.polanski@warwick.ac.uk>

# this is where we start
ENTRYPOINT ["bash", "/scripts/demo/gp2s_tarwrapper.sh"]