#note the test.bed - made out of the .broadPeak file by just taking the first three columns
docker run -it --rm -v /home/kpolanski/docker_demo/testdata-k562:/agave wellington test.bed wgEncodeUwDgfK562Aln.bam
