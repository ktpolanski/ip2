# run as:
#   docker build -f csi/Dockerfile  .

# base everything on a recent Ubuntu
FROM debian:latest

# get system packages up to date then install a basic scientific python
RUN apt-get update && apt-get -y upgrade && \
    apt-get -y install python3 \
         python3-numpy python3-scipy python3-pandas python3-matplotlib \
	 python3-seaborn python3-h5py

# add and configure our code
RUN mkdir /scripts
COPY csi/ html/ /scripts/
MAINTAINER Sam Mason <sam@samason.uk>
ENTRYPOINT ["bash", "/scripts/csi/csi_wrapper.sh"]
CMD ["--help"]
