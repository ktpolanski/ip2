docker run -it --rm -v /home/kpolanski/docker_demo/TAIR10:/agave -w /agave hmt-index genome_rm.fa annot.gff3 gene_id 500 minimeme.txt 0.05 5 --No --No

docker run -it --rm -v /home/kpolanski/docker_demo/TAIR10:/agave -w /agave hmt-index genome_rm.fa annot.gff3 gene_id 500 MeMe_friendly_joint.txt 0.05 5 --No --No

docker run -it --rm -v /home/kpolanski/docker_demo/TAIR10_nooverlap:/agave -w /agave hmt-index genome_rm.fa annot.gff3 gene_id 500 MeMe_friendly_joint.txt 0.05 5 --No --NoOverlap