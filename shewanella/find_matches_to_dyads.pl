#!/usr/bin/perl -w

#############################################################################
#                                                                           #
#   copyright            : (C) 2000-2003 by David J. Studholme              #
#   email                : ds2@sanger.ac.uk                                 #
#                                                                           #
#############################################################################

#############################################################################
#*                                                                          #
#*   This program is free software; you can redistribute it and/or modify   #
#*   it under the terms of the GNU General Public License as published by   #
#*   the Free Software Foundation; either version 2 of the License, or      #
#*   (at your option) any later version.                                    #
#*                                                                          #
#############################################################################

use warnings;
use strict ;
use Bio::SeqIO;


my $genome_sequence_file = shift or die "Usage: $0 <genome sequence file (.fna)>  <ptt file>\n" ;
my $ptt_file = shift or die "Usage: $0 <genome sequence file (.fna)>  <ptt file>\n" ;

my $dirname = "." ;

### Read the dyad sequences
my %dyads;
opendir(DIR, "$dirname") or die "can't opendir $dirname: $!";
while (defined(my $dyad_file = readdir(DIR))) {
    if ( $dyad_file =~ m/significant-dyads.fna$/ ) {
	warn "Reading sequence from file '$dyad_file'\n";
	my $inseq = Bio::SeqIO->new('-file' => "<$dyad_file",
				    '-format' => 'fasta' ) ;
	while (my $seq_obj = $inseq->next_seq ) {
	    my $id = $seq_obj->id;
	    my $sequence = $seq_obj->seq;
	    my $desc = $seq_obj->description;
	    $dyads{$sequence} = $id;
	}
    }
} 
closedir(DIR);

warn "Got ".(keys %dyads)." dyad sequences\n";

### Read the genome sequence
warn "Reading sequence from file '$genome_sequence_file'\n";
my ($sequence, $id, $desc);
my $inseq = Bio::SeqIO->new('-file' => "<$genome_sequence_file",
			    '-format' => 'fasta' ) ;
while (my $seq_obj = $inseq->next_seq ) {
    $id = $seq_obj->id;
    $sequence = $seq_obj->seq;
    $desc = $seq_obj->description;
}
my $seqlen = length($sequence) ;
warn  "$id ($desc) is $seqlen residues long\n\n" ; 

### Get the intergenic sequences
my @intergenic_seqs = get_intergenic_sequences($sequence, $ptt_file);
warn "Got ".@intergenic_seqs." intergenic sequences\n";

### Find all the matches to all the dyads
my %dyad2matches;
foreach my $dyad (sort keys %dyads) {
    foreach my $intergenic_seq (@intergenic_seqs) {
	#warn "$intergenic_seq =~ m/($dyad)/\n";
	if ($intergenic_seq =~ m/($dyad)/) { 
	    #warn "$dyad matches $1\n";
	    $dyad2matches{$dyad}{$1} = 1;
	}
	
	if (revcomp($intergenic_seq) =~ m/($dyad)/) {	
	    #warn "$dyad matches $1\n";
	    $dyad2matches{$dyad}{$1} = 1;
	}
    }
}


### Do any of the dyads overlap in their matching sequences?
foreach my $dyad1 (sort keys %dyad2matches) {
    foreach my $dyad2 (sort keys %dyad2matches) {
	unless ($dyad1 eq $dyad2) {
	    foreach my $seq1 (keys %{$dyad2matches{$dyad1}} ){
		foreach my $seq2 (keys %{$dyad2matches{$dyad2}} ){
		    if ($seq1 =~ m/$seq2/ or
			$seq2 =~ m/$seq1/ ) { 
			#warn "$seq1 matches $dyad1 and $dyad2\n";
			if ( (keys %{$dyad2matches{$dyad1}}) > (keys %{$dyad2matches{$dyad2}}) ) {
			    delete $dyad2matches{$dyad2};
			} else {
			    delete $dyad2matches{$dyad1};
			}
		    }
		}
	    }
	}
	
    }
}


foreach my $dyad (sort keys %dyad2matches) {

    my $clustal_path = "/home/studhold/clustalw1.83/";
    my $hmmer_path = "/home/studhold/hmmer-1.8.5/";
    
    my $n = keys %{$dyad2matches{$dyad}};
    if ($n > 2) {
    
	### make an alignement from each set of dyad-matches
	my $tmp_seq_file = "$dyad.matches.fna";
	open(TMP,">$tmp_seq_file") or die "Failed to open file for writing '$tmp_seq_file' $!\n";
	
	my $i=0;
	foreach my $seq (sort keys %{$dyad2matches{$dyad}}) {
	    $i++;	
	    print TMP ">$i\n$seq\n";
	}
	close TMP;

	my $align_file = "$dyad.matches.msf";
	my $hmmsearch_file = "$genome_sequence_file.$align_file.hmmsearch";
	my $hmm_file =  "$align_file.hmm";
	
	unless (-s $align_file and
		-s $hmmsearch_file and
		-s $hmm_file) {
	    
	    
	    my $align_cmd = "$clustal_path/clustalw $tmp_seq_file -OUTPUT=GCG -endgaps -gapopen=100";
	    my $hmm_build_cmd = "$hmmer_path/hmmb $hmm_file $align_file";
	    my $hmm_search_cmd = "$hmmer_path/hmmls -t 10 -c -F $hmm_file $genome_sequence_file > $hmmsearch_file";
	    
	    my $cmd = "$align_cmd && $hmm_build_cmd && $hmm_search_cmd";
	    my $lsf_cmd = "bsub -q short -o lsf.out '$cmd'";
	    
	    
	    if (0) {
		warn "$cmd\n\n";
		system($cmd);
		
	    } else {
		warn "$lsf_cmd\n\n";
		system($lsf_cmd);
		
	    }
	}
	
	
	
    }
}



exit ;


sub read_ptt_line {
  # Expects a scalar value containing one line of a .ptt file
  # Extracts the 'location', strand, and 'product'
  #   from that line and returns them, in that order
  
  my $read_line = $_[0];
  my $i;
  my $product;
  my @fields;
  
  if ($read_line =~ /(\d+)\.\.(\d+)\s+([+-])/){
    # Extract 'location' and 'strand'    
    @fields = split(/\s+/,$read_line);
    $i = 5;
    $product = "";
    while ($fields[$i]) {
      # Extract 'product'
      $product = $product." ".$fields[$i];
      $i++;
    } 
    return ($1, $2, $3, $product);
  } else {
    return undef ;
  }
  
} ; # end of sub read_ptt_line;


sub get_intergenic_sequences {
  my $minimum_size = 50 ;
  my $maximum_size = 300;
  my $sequence = shift or die;
  my $ptt_file = shift or die;
  my @intergenic_sequences;


  open(PTT_FILE, $ptt_file) or  die "Could not open .ptt file: $ptt_file\n" ;
  print STDERR "\nOpened the .ptt file\n" ;
  
  my ($start_pos1, $end_pos1, $strand1, $product1) = ( 0, 0, "+", "Origin" ) ;
  while (my $read_line = <PTT_FILE>) {  
    my @fields = read_ptt_line($read_line) ;
    if ( @fields > 1  ) {
   
      my ($start_pos2, $end_pos2, $strand2, $product2) = @fields ;
      
      my $intergenic_sequence ;
      
      my $size = $start_pos2 - $end_pos1  ;
      if ( $size >= $minimum_size ) {
	
	### Ignore region if it is only at 3' end of genes
	if ( $strand1 eq "+" and $strand2 eq "-" ) {
	  $intergenic_sequence = "" ;
	} 
	
	### If intergenic region is small, then no problem
	elsif ( $size <= $maximum_size ) {
	  $intergenic_sequence = substr($sequence, $end_pos1, $size) ;
	}	
	
	### If both genes are on +ve strand
	elsif ( $strand1 eq "+" and $strand2 eq "+" ) { 
	  $intergenic_sequence = substr($sequence, $start_pos2-$maximum_size, $maximum_size) ;	
	}
	
	### If both genes are on +ve strand
	elsif ( $strand1 eq "-" and $strand2 eq "-" ) { 
	  $intergenic_sequence = substr($sequence, $end_pos1, $maximum_size) ;	
	}

	### If genes are divergent
	elsif ( $strand1 eq "-" and $strand2 eq "+" ) { 
	  my $half_max = $maximum_size / 2 ;
	  my $intergenic_sequence1 = substr($sequence, $end_pos1, $half_max) ;	
	  my $intergenic_sequence2 = substr($sequence, $start_pos2-$maximum_size, $half_max) ;	
	  $intergenic_sequence = $intergenic_sequence1.$intergenic_sequence2 ;
	}
	
	else {
	  die "What the hell happened here then?\n" ;
	}
	
	### OK, add the intergenic sequence to the list!
	push @intergenic_sequences, $intergenic_sequence  ;
	
      } 
      
      #my $length = length($intergenic_sequence) ;
      #print "length = $length\n" ;

      ($start_pos1, $end_pos1, $strand1, $product1) = ($start_pos2, $end_pos2, $strand2, $product2)
    }
    else {
      #print "Ignoring line: $read_line\n" ;
    }

  }
  return @intergenic_sequences ;
}



sub revcomp {
  my $string = shift @_ ;
  my $seq = "" ;
  my @array = split(//, $string) ;
  while (my $char = pop @array) {
    $char =~ tr/ACGT/TGCA/ ;
    $seq = $seq.$char 
  }
  #print STDERR "Revcomp of $string -> $seq\n" ;
  return $seq 
}
