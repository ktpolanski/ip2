#(once inside the bogaotory/docker-perl-bio image)
source /root/.profile
perl get_significant_dyads.pl shewanella_perldyads.fasta annotation.ptt > get_significant_dyads_out.txt
cpanm Graph::Undirected
