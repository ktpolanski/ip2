FROM ubuntu:latest

RUN apt-get update && apt-get -y install software-properties-common && add-apt-repository multiverse
RUN apt-get update && apt-get -y upgrade
RUN apt-get -y install build-essential g++ libboost-all-dev wget cmake git pkg-config nvidia-cuda-toolkit
RUN wget http://bitbucket.org/eigen/eigen/get/3.2.5.tar.gz && tar -xf 3.2.5.tar.gz
RUN mkdir eigenbuild && cd eigenbuild && cmake /eigen-eigen-bdd17ee3b1b3 && make install
RUN cd / && git clone https://github.com/smason/mdipp && cd mdipp && make ncuda=1

ENTRYPOINT ["/mdipp/mdi++"]