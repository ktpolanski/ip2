# run as:
#   docker build -t hmt .
# (with hmt.py and hmt_wrapper.sh in the current directory)

# base everything on a recent Ubuntu
FROM debian:latest

# get system packages up to date then install a basic scientific python
RUN apt-get update && apt-get -y upgrade && \
    apt-get -y install python3 \
         python3-numpy python3-scipy python3-pandas \
         ttf-bitstream-vera

# add the code
RUN mkdir /scripts
COPY hmt.py hmt_wrapper.sh /scripts/

MAINTAINER Krzysztof Polanski <k.t.polanski@warwick.ac.uk>

# so this is what is going to run by default when you trigger this, in the virtual machine
# call the wigwams wrapper from the other directory while staying in /agave with the files
ENTRYPOINT ["bash", "/scripts/hmt_wrapper.sh"]