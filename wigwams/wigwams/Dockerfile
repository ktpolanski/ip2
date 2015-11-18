# run as:
#   docker build -t wigwams .
# (with wigwams.py and wigwams_wrapper.py in the current directory)

# base everything on a recent Ubuntu
FROM debian:latest

# get system packages up to date then install a basic scientific python
RUN apt-get update && apt-get -y upgrade && \
    apt-get -y install python3 \
         python3-numpy python3-scipy python3-pandas python3-matplotlib python3-seaborn \
         ttf-bitstream-vera

# add and configure our code

# take wigwams.py and wigwams_wrapper.py, which we need
# and put it all into /wigwams inside the docker container
RUN mkdir /wigwams
COPY wigwams.py wigwams_wrapper.py /wigwams/

MAINTAINER Krzysztof Polanski <k.t.polanski@warwick.ac.uk>

# so this is what is going to run by default when you trigger this, in the virtual machine
# call the wigwams wrapper from the other directory while staying in /agave with the files
ENTRYPOINT ["python3", "/wigwams/wigwams_wrapper.py"]

# if nothing else is specified in the docker call, just run --help
CMD ["--help"]