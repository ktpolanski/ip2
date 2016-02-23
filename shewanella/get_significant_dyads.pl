#!/usr/bin/perl

### Find statistically significant dyads in intergenic regions

#############################################################################
#                                                                           #
#   copyright            : (C) 2007 by David J. Studholme                   #
#   email                : david.studholme@tsl.ac.uk                        #
#                                                                           #
#############################################################################

#############################################################################
#                                                                           #
#    This program is free software; you can redistribute it and/or modify   #
#    it under the terms of the GNU General Public License as published by   #
#    the Free Software Foundation; either version 2 of the License, or      #
#    (at your option) any later version.                                    #
#                                                                           #
#############################################################################

# Inspired by: 
# van Helden J, Rios AF, Collado-Vides J. 
# Discovering regulatory elements in non-coding sequences by analysis of spaced dyads.
# Nucleic Acids Res. 2000 Apr 15;28(8):1808-18.
# PMID: 10734201 
# and 
# Li H, Rhodius V, Gross C, Siggia ED. 
# Identification of the binding sites of regulatory proteins in bacterial genomes.
# Proc Natl Acad Sci U S A. 2002 Sep 3;99(18):11772-7.
# PMID: 12181488

use strict ;
use POSIX qw(log10);
use Math::BigInt ;
use Math::BigFloat ;
use Bio::SeqIO ;

### Word length
my $k = 4 ;
warn "Word length, k = $k\n" ;

### Spacer lengths
my $minspace = 10 ;
my $maxspace = 20 ;
warn "Space lengths $minspace - $maxspace\n" ;

### Genome files
my $genome_sequence_file = shift or die "Usage: $0 <seq_file> <ptt_file>\n" ;
my $ptt_file = shift or die "You must specify a ptt file\n";

### Read the seqfile into a string
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

### Pre-compute values of factorials to save compute time later
warn "Pre-computing factorials\n";
my $precalc = precalculate_factorials(4000) ;
warn "done\n\n" ;

### Now get a list of intergenic sequences
warn "Getting intergenic sequences\n";
my @intergenic_sequences = get_intergenic_sequences($sequence, $ptt_file) ;
warn "Finished getting intergenic sequences\n";

### Get a list of all possible words of length $k
my @words = generate_all_words($k, $k) ;

### Np is the number of possible dyads considered in each analysis
my $number_of_words = @words ;
my $Np = $number_of_words * $number_of_words * (1 + $maxspace - $minspace) ;
warn "Number of possible words = $number_of_words\n" ;
warn "Number of possible dyads, Np  = $Np\n\n" ;

### Stringency of statistical significance
#my $threshold = 0.001 / $Np;
my $threshold = 0.000000001;
warn "Significance threshold = $threshold \n" ;

### Find number of non-overlapping occurrences of each word
warn "Calculating word frequencies in background sequence\n" ;
my %word_freqs  ; #  hash of observed frequencies for each word
foreach my $word (@words) {
  my $freq =  word_freq($word, @intergenic_sequences) ;
  #warn "$word: $freq \n" ;
  $word_freqs{$word} = $freq ;
} 
warn "\nFinished calculating word frequencies in background sequence\n\n" ;

my %significant_dyads ;

# ref to a hash of observed frequencies for each dyad
my $dyad_occs = {} ; 

# record number of positions that were not counted because of overlaps
my $excluded_positions = {} ;

# The number of possible positions is dependent on the dyad length
my %T_count  ; 

### For each spacer size and for (all the intergenic) sequences calculate dyad frequencies ;
warn "Calculating frequencies of each dyad.\n";
for (my $spacer = $minspace; $spacer <= $maxspace ; $spacer++) {
  
  # record number of positions that were not counted because of overlaps
  my $excluded_positions = {};

  ### First scan the forward strands
  foreach my $intergenic_sequence (@intergenic_sequences) { 
    ($dyad_occs, $excluded_positions, $T_count{$spacer})  = 
      count_dyads($dyad_occs, $excluded_positions, $intergenic_sequence, $T_count{$spacer}, $k, $spacer) ;
  }
  
  ### Now scan the reverse strands
  foreach my $intergenic_sequence (@intergenic_sequences) { 
    ($dyad_occs, $excluded_positions, $T_count{$spacer})  = 
      count_dyads($dyad_occs, $excluded_positions, revcomp($intergenic_sequence), $T_count{$spacer}, $k, $spacer) ;
  }
}
warn "Finished calculating frequencies of each dyad.\n\n";

### examine each dyad for over-rpresentation
foreach my $dyad (sort keys %$dyad_occs) {
    
    ##warn "Examining dyad '$dyad'\n";
    
    ### Some positions in the sequence are excluded to avoid counting overlapping occurrences of dyads
    $$excluded_positions{$dyad} = 0 unless $$excluded_positions{$dyad} ;
    
    ### Parse the dyad sequence
    $dyad =~ m/([ACGT]{$k})(.*)([ACGT]{$k})/i ;
    my ($word1, $word2) = ($1, $3) ;
    my $spacer = length($2) ;
    
    ### Expected frequency of dyad is product of probabilities for each of the constituent words
    my $expected_dyad_freq = $word_freqs{$word1} * $word_freqs{$word2} ;
    
    ### Calculate the actual number of possible positions, taking into account any excluded positions
    my $T_adjusted = $T_count{$spacer} - $$excluded_positions{$spacer} ;
    
    ### Calculate the expected number of occurrences of this dyad in the sequence
    my $expected_dyad_occ = $expected_dyad_freq * $T_adjusted ;
    
    ### Calculate P value - a measure of statistical significance for this dyad
    my $P = 1 ;
    $P = P_value($T_adjusted, $word_freqs{$word1}, $word_freqs{$word2}, $precalc, $$dyad_occs{$dyad} ) if $$dyad_occs{$dyad} > $expected_dyad_occ ;
    
    
    
    
    ### is this dyad is statistically significant?
    if ( $P <= $threshold ) {
	my $rdyad = revcomp($dyad) ;
	warn "SIGNIFICANT:\t$dyad / $rdyad:\t$$dyad_occs{$dyad} occurrences, $expected_dyad_occ expected, P(X>=$$dyad_occs{$dyad}) = $P\n";
	$significant_dyads{$dyad} = $P ;
    } else {
	#my $rdyad = revcomp($dyad) ;
	#warn "NOT SIGNIFICANT:\t$dyad / $rdyad:\t$$dyad_occs{$dyad} occurrences, $expected_dyad_occ expected, P(X>=$$dyad_occs{$dyad}) = $P\n";
    }
     
} 

my $dyad_count = (keys %significant_dyads) ;
warn "\nNumber of statistically significant dyads: $dyad_count\n\n" ;


 ### print the significant dyads
my $i = 1;
foreach my $dyad_seq (sort keys %significant_dyads) {
    my $P = $significant_dyads{$dyad_seq};
    print ">dyad_$i P=$P\n$dyad_seq\n";
    $i++
}

exit ;

  


sub get_intergenic_sequences {
  my $minimum_size = 50 ;
  my $maximum_size = 300 ;
  my ($sequences, $ptt_file) =  @_ ;
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
      
      my $length = length($intergenic_sequence) ;
      #print "length = $length\n" ;

      ($start_pos1, $end_pos1, $strand1, $product1) = ($start_pos2, $end_pos2, $strand2, $product2)
    }
  }
  return @intergenic_sequences ;
}





















sub word_freq {
  ### Return the frequency of non-overlapping occurrences of $word in @sequences
  ### Scan both strands

  my ($word, @sequences) = @_ ;
  my $word_count = 0 ;
  my $T_count = 0 ; # number of possible positions
  my $k = length($word) ;
  
  ### Scan forward strand
  foreach my $sequence (@sequences) {
    my $T = length($sequence) ;
    for (my $i = 0 ; $i <= ($T-$k) ; $i++ ) {
      $T_count++; # number of possible positions
      if ( $word eq substr($sequence, $i, $k) ) {
	$word_count++ ; # number of occurrences
	$i = $i + ($k-1) ; # move along to avoid overlapping occurrences
      }
    }
  }

  ### Scan reverse strand
  foreach my $sequence (@sequences) {
    my $T = length(revcomp($sequence)) ;
    for (my $i = 0 ; $i <= ($T-$k) ; $i++ ) {
      $T_count++; # number of possible positions
      if ( $word eq substr($sequence, $i, $k) ) {
	$word_count++ ; # number of occurrences
	$i = $i + ($k-1) ; # move along to avoid overlapping occurrences
      }
    }
  }

  my $frequency = $word_count / $T_count ;
  return $frequency ;
}


sub count_dyads {
  my ($dyad_occs, $excluded_positions, $sequence, $T_count, $k, $spacer) = @_ ; 
  my $T = length($sequence) ;
  
  my $window_size = $k + $spacer + $k ;

  my %previous_match_to_this_dyad ; # to keep track of overlapping matches

  for ( my $i = 0 ; $i <= ($T-$window_size) ; $i++ ) {
    $T_count++ ; # number of positions checked for presence of dyad
    my $word1 = substr($sequence, $i, $k) ;
    my $word2 = substr($sequence, $i+$k+$spacer, $k) ;
    my $dyad = $word1 ;
    for (my $i = 1; $i <= $spacer ; $i++) {
      $dyad .= "." ;
    }
    $dyad .= $word2;

    $previous_match_to_this_dyad{$dyad} = (-1  - $k)  unless ($previous_match_to_this_dyad{$dyad}) ;
    if ( $previous_match_to_this_dyad{$dyad} >= ($i-$k) ) {
      ### The previous occurrence of this dyad overlaps this occurrence
      $$excluded_positions{$dyad}++ ;
    } else {
      ### OK - no overlap with previous occurrence
      $$dyad_occs{$dyad}++ ;
      $previous_match_to_this_dyad{$dyad} = $i ;
    }
  }
  return ($dyad_occs, $excluded_positions, $T_count) ;
}


sub generate_all_words {
  ### Return every possible word

  my ($minword, $maxword) =  @_ ;
  my @words = ( "", "", "", "" ) ;
  my %words ;

  for (my $word_length = $minword; $word_length <= $maxword; $word_length++) { 
    my $count = 0 ;
    until ( $count == $word_length ) {
      @words = extend_words(@words) ;
      $count++ ;
    }
  }

  ### Make the list non-redundant
  foreach my $word (@words) {
    $words{$word} = 1 
  }
  
  return (sort keys %words)
}


sub extend_words {
  my @old_words = @_ ;
  my @new_words ;
  foreach my $letter ("A", "C", "G", "T" ) {
    foreach my $word ( @old_words ) {
      push @new_words, "$word$letter" 
    }
  }
  return ( @new_words )
}


  
sub P_value {
  ### Calculate P = P(X >= N) where N = 0, 1, 2, 3, ...
  ###  and assuming a Poisson distribution
  ### Note that P = 1 - P_ where P_ = P(X<N)
  my ($T, $fW1, $fW2, $precalc, $N)  = @_ ;
  my $mu = $T * $fW1 * $fW2 ;
  my $P_ = 0 ;
  my $exp = Math::BigFloat->new("1") ; 
  $exp *= exp(-$mu) ; 
  for (my $i = 0 ; $i < $N ; $i++) {
    my $tmp = ($mu ** $i) * $exp / factorial($i, $precalc) ; 
    $tmp = 0 unless $tmp =~ /[d+\.]+/ ;
    #print STDERR "mu=$mu, P(X=$i) = $tmp\n" ;
    $P_ += $tmp ;
  }
  my $P = 1 - $P_ ;
  $P =  sprintf("%.15f", $P) ;
  #print STDERR "\nP(X>=$N) = $P\n" ;
  return $P ;
}



sub factorial {
  my ($n, $precalc) = @_ ;
  return 1 if $n < 2 ;
  my $result = $$precalc{$n}  ;
  die "No precomputed value for $n factorial\n" unless $result ;
  return $result ;
}


sub precalculate_factorials {
  ### Pre-calculate factorials, which are computationally expensive
  my $max = shift @_ ;
  my $result = Math::BigFloat->new("1") ;
  my %factorial ;
  for (my $i = 1 ; $i <= $max ; $i++) {
    $result *= $i ;
    $factorial{$i} = $result ;
    #print STDERR "factorial($i) = $result\n" ;
  } 
  return \%factorial
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



sub get_score {
  ### Find the maximum alignment score between two dyads
  my ($dyad1, $dyad2)  = @_ ;
  my ($long_dyad, $short_dyad) ;
  my $max_score = 0 ;

  ### Which is the longer dyad?
  if ( length($dyad1) > length($dyad2) ) {
    $long_dyad = $dyad1 ;
    $short_dyad = $dyad2 ;
  } else {
    $long_dyad = $dyad2 ;
    $short_dyad = $dyad1 ;
  }
  
  ### The short dyad can be aligned against the long dyad
  ###  at any position with an offset of between (1-short) and (long-1)
  for (my $i = (1 - length($short_dyad)) ; $i <= (length($long_dyad) + 1) ; $i++ ) {
    
    ### Calculate the score
    my $score = 0 ;
    for (my $j = 0 ; $j < length($short_dyad) ; $j++) {
      if ( ($j+$i > -1) and (($j+$i) < length($long_dyad)) ) {
	my $char1 = substr($short_dyad, $j, 1) ;
	my $char2 = substr($long_dyad, $j+$i, 1) ;

	if ($char1 =~ m/[ACGT]/ and $char2 =~m/[ACGT]/) { 
	  $score++ if ($char1 eq $char2) ;
	  $score-- if ($char1 ne $char2) ;
	}
      }
    }
    ### Record this score if it is the highest
    $max_score = $score if $score > $max_score ;
  }
  #print "best score for $dyad1 vs $dyad2: $max_score\n" ;
  return $max_score ;
}

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

