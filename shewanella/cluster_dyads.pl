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

### Classify dyads into clusters


# Inspired by: 

# Li H, Rhodius V, Gross C, Siggia ED. 
# Identification of the binding sites of regulatory proteins in bacterial genomes.
# Proc Natl Acad Sci U S A. 2002 Sep 3;99(18):11772-7.
# PMID: 12181488



use strict ;
use Graph::Undirected ;

### clustering theshold
my $build_threshold = 7 ;
my $maintain_threshold = 5 ;
print "Clustering thresholds: $build_threshold, $maintain_threshold\n" ;

### Read all the significant dyads from file
my %significant_dyads ;
my $infile = shift or die "Usage: $0 <file>\n" ;
open(FILE, "<$infile") or die "Could not open dyad infile: $infile\n" ;
while ( my $readline = <FILE> ) {
  if ( $readline =~ m/([ACGTN]+)\s\([ACGTN]+\)/ ) {
    my $dyad = $1 ;
    print STDERR "Read: $dyad\n" ;
    $significant_dyads{$dyad} = 1 ;
  }
}

my $dyad_count = (keys %significant_dyads) ;
print "\nNumber of statistically significant dyads: $dyad_count\n\n" ;


### Generate an undirected Graph object where vertices represent dyads
###  and edges connect similar ones
print STDERR "Creating graph object ... " ;
my $graph = Graph::Undirected->new(keys %significant_dyads) ;
print STDERR "done\n\n" ;


### Compare every dyad against every dyad
my ($total_score, $score_count) = (0, 0) ;
print STDERR "\nComparing all against all dyads ... " ;
my @dyads = (sort keys %significant_dyads) ;
while ( my $dyad2 = shift @dyads ) {
  foreach my $dyad1 ( @dyads ) {
    my $score = get_score($dyad1, $dyad2) ;
    $total_score +=  $score;
    $score_count++;
    if ( $score >= $build_threshold ) {
      $graph->add_edge($dyad1, $dyad2)   ;
      #print STDERR "Added edge $dyad1 - $dyad2 (score = $score )\n"  ;
    }
  }
}
print STDERR "done\n\n" ;


### Calculate average dyad similarity
my $average = $total_score / $score_count ;
print "Average similarity score between dyads = $average\n" ;
print STDERR "Average similarity score between dyads = $average\n" ;

$graph = refine_clusters($graph) ;

print_graph($graph, "graphfile") ;

### Derive final clusters of dyads from the graph object
my @clusters = $graph->strongly_connected_components() ;
my $number_of_clusters = @clusters ;
print "Final number of clusters: $number_of_clusters\n" ;

### Print out these clusters to fasta files 
my $cluster_num = 1 ;
foreach my $cluster (@clusters) {
  my @vertices = @$cluster ;
  my $number_of_vertices = @vertices ;
  if ($number_of_vertices > 0 ) {
    open(OUTFILE, ">dyads_cluster\_$cluster_num.fa") or die "failed to open outfile\n" ;
    my $dyad_number = 1 ;
    foreach my $vertex (@vertices) {
      my $edges = $graph->edges($vertex) ;
      print OUTFILE ">$cluster_num\_$dyad_number\n$vertex\n" ;
      $dyad_number++ ;
    }
    close OUTFILE ;
  }
  $cluster_num++ ; 
}


exit ;


sub print_graph {
  my ($graph, $outfile) = @_ ;

  open (GRAPH, ">$outfile") or die "Failed to open $outfile\n" ;

  ### Compare every dyad against every dyad
  my ($total_score, $score_count) = (0, 0) ;
  print STDERR "\nCreating graphfile ... " ;
  my @dyads = (sort keys %significant_dyads) ;
  while ( my $dyad2 = shift @dyads ) {
    print "$dyad2\n" ;
    foreach my $dyad1 ( @dyads ) {
      my $score = get_score($dyad1, $dyad2) ;
      $total_score +=  $score;
      $score_count++;
      if ( $graph->has_edge($dyad1, $dyad2) ) {
	print GRAPH "$dyad1 $dyad2\n"  ;
      } 
    }
  }
  close GRAPH ;
  print STDERR "done\n\n" ;
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
  
  ### Convert the dyad sequences to arrays
  my @short_dyad = split("", $short_dyad) ;
  my @long_dyad = split("", $long_dyad) ;
  
  ### The short dyad can be aligned against the long dyad
  ###  at any position with an offset of between -4 and (long + 4 - short)
  for (my $i = -4 ; $i <= (length($long_dyad) + 4 - length($short_dyad)) ; $i++ ) {
    
    ### Calculate the alignment score for this offset
    my $score = 0 ;
    for (my $j = 0 ; $j < length($short_dyad) ; $j++) {
      
      if ( ($j+$i > -1) and 
	   (($j+$i) < length($long_dyad)) and
	   ($j < length($short_dyad)) ) {
	
	my $char1 = $short_dyad[$j] ;
	my $char2 = $long_dyad[$j+$i] ;
	
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


sub refine_clusters {

  my ($graph) = @_ ;
  
  print STDERR "Refining (breaking up) large clusters " ;
    
  my @clusters = $graph->strongly_connected_components() ;
  
  my $number_of_clusters = @clusters ;
  print STDERR "\n\nnumber of clusters: $number_of_clusters\n" ;
  
  my $cluster_num = 1 ;
  
  foreach my $cluster (@clusters) {
    #print STDERR "\n\n### Cluster $cluster_num ##\n" ;
    my @vertices = @$cluster ;
    my $number_of_vertices = @vertices ;
    my %cluster_members ;
    foreach my $member (@vertices) {
	$cluster_members{$member} = 1 ;    
      }


    ### Compare every dyad against every dyad within this cluster
    my ($total_score, $score_count) = (0, 0) ;
    print STDERR "Comparing dyads within this cluster $cluster_num ... " ;
    while ( my $dyad2 = shift @vertices ) {
      foreach my $dyad1 ( @vertices ) {
	my $score = get_score($dyad1, $dyad2) ;
	print STDERR "  $dyad1 : $dyad2 score = $score\n" ;
	print STDERR "  Cluster should be severed somewhere between $dyad1 and $dyad2\n" if $score < $maintain_threshold ;
	my @dyad1_neighbours = $graph->neighbours($dyad1) ;
	my @dyad2_neighbours = $graph->neighbours($dyad2) ;
	if ( @dyad1_neighbours < @dyad2_neighbours ) {
	  foreach my $neighbour (@dyad1_neighbours) {
	    $graph->delete_edge($dyad1, $neighbour) if $cluster_members{$neighbour} ;
	    print STDERR "Deleted edge $dyad1 - $neighbour\n"  if $cluster_members{$neighbour} ;;
	  } 
	} else {
	  foreach my $neighbour (@dyad2_neighbours) {
	    $graph->delete_edge($dyad2, $neighbour)  if $cluster_members{$neighbour} ;;
	    print STDERR "Deleted edge $dyad2 - $neighbour\n"  if $cluster_members{$neighbour} ;;
	  } 
	}
	
      }
    }
    close GRAPH ;
    print STDERR "done\n" ;
    
    $cluster_num++ ; 
  }
  print STDERR "done.\n\n" ;
  return $graph ;
}
