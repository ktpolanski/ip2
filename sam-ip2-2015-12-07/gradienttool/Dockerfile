# run as:
#   docker build -f Dockerfile -t gradienttool .

# base everything on a recent Ubuntu
FROM debian:latest

# get system packages up to date then install a basic python
# note the prerequisites for scipy updating - liblapack-dev and gfortran
# in turn, matplotlib when pipped needs pkg-config and libfreetype6-dev
RUN apt-get update && apt-get -y upgrade && \
    apt-get -y install python3 python3-dev python3-setuptools liblapack-dev gfortran build-essential \
    pkg-config libfreetype6-dev ttf-bitstream-vera

# set up pip and grab stuff through that, to make all the installations homogeneous
RUN easy_install3 pip
RUN pip install cython numpy scipy matplotlib pandas seaborn
# this needs to be separate as the primadonna of a package needs to see numpy installed
RUN pip install GPy

# set up code
RUN mkdir /scripts
COPY scripts/ /scripts/
ENTRYPOINT ["bash", "/scripts/gradienttool_tarwrapper.sh"]