#note the K562 data has been moved to the bootstrap folder, if you need to re-run the footprint for some reason
docker run -it --rm -v /home/kpolanski/docker_demo/testdata-k562:/agave wellington-footprint wgEncodeUwDgfK562Hotspots.broadPeak wgEncodeUwDgfK562Aln.bam

docker run -it --rm -v /home/kpolanski/docker_demo/testdata-boot:/agave wellington-bootstrap wgEncodeUwDgfHepg2Hotspots.broadPeak wgEncodeUwDgfK562Hotspots.broadPeak --start wgEncodeUwDgfHepg2Aln.bam wgEncodeUwDgfK562Aln.bam

docker run -it --rm -v /home/kpolanski/docker_demo/testdata-bootmini:/agave wellington-bootstrap ex_reg.bed --start ex1.bam ex2.bam