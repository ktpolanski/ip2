# run as:
#   docker build -f Dockerfile-gradient -t gradienttool .

# base everything on a recent Ubuntu
FROM debian:latest

# get system packages up to date then install a basic scientific python
# note the prerequisites for scipy updating - liblapack-dev and gfortran
RUN apt-get update && apt-get -y upgrade && \
    apt-get -y install python3 python3-dev python3-setuptools liblapack-dev gfortran build-essential \
    python3-numpy python3-scipy python3-matplotlib

# set up pip and grab GPy and pandas through that
RUN easy_install3 pip
RUN pip install GPy pandas

# set up code
RUN mkdir /gradienttool
COPY gradienttool/ /gradienttool/
ENTRYPOINT ["python3", "/gradienttool/__main__.py"]
CMD ["--help"]