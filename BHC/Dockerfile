#run as:
#	docker build -t bhc .

#grab existing R image
FROM r-base:latest

#need Python for R's argparse port
RUN apt-get -y update && apt-get -y upgrade
RUN apt-get -y install python3 ttf-bitstream-vera

#copy over scripts
RUN mkdir /scripts
COPY scripts /scripts

#R-side setup
RUN Rscript /scripts/setup.R

#that's it, ready to use
ENTRYPOINT ["Rscript", "/scripts/bhc.R"]