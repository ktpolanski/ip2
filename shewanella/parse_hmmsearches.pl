#!/usr/bin/perl

use strict ;
use warnings ;
use Bio::SearchIO;
use Bio::SeqIO;
use CGI;
use Graph::Directed;
use Cwd;

my $operon_threshold = 75;
my $default_threshold = 18;
my $cgi = new CGI;
my $significance_threshold = 0.001;
my $offset = 200; # allow intergenic region to extend into the annotated gene due to poor start codon prediction


my $usage = "Usage: $0 <PTT file>";
my $ptt_file = shift or die "$usage\n";
my $dirname = ".";

my $genome;
if ($ptt_file =~ m/^(\S+)\.ptt$/i) {
    $genome = $1;
} else {
die "$usage\n"
}



### Read in all the clusters into hash
my %dyad2cluster ;

my @hmmsearch_files;
opendir(DIR, "$dirname") or die "can't opendir $dirname: $!";
while (defined(my $file = readdir(DIR))) {
    
    if ( $file =~ m/^$genome(\S+)\.hmmsearch$/) {
	push @hmmsearch_files, $file;
    }
}


die "There are too many (".@hmmsearch_files.") hmmsearch output files.\n" if @hmmsearch_files > 500;





print $cgi->start_html(-title=>"Search results for $genome",
		       -author=>'david.studholme@sainsbury-laboratory.ac.uk',
		       -meta => {'keywords'=>'',
				 'copyright'=>'David Studholme'},
		       -text => 'black',
		       -bgcolor => 'white',
		       ); 	

print "<link href=\"http://www.ayeaye.tsl.ac.uk/templates/sainsburylaboratory/css/template_css.css\" rel=\"stylesheet\" type=\"text/css\"/>\n";


my $cwd = getcwd();
print "<h1>$cwd</h1>\n";

### A graph object to store predicted operon structure
my $operons = new Graph::Directed;
	
### Read the .PTT file
my ($start_1, $end_1, $strand_1, $id_1, $desc_1) = (0,0,"", "", "");
my (%upstream_start, %upstream_end, %id2desc, %id2strand, %gene_start, %gene_end);
open (FILE, "<$ptt_file") or die "Failed to open $ptt_file\n$!";
while (<FILE>) {
    chomp;
    if ( m/^
	 (\d+)
	 \.\.
	 (\d+)
	 \t
	 ([+-])
	 \t
	 \S+
	 \t
	 \S+
	 \t	 
	 (\S+)
	 \t
	 (\S+)
	 \t
	 \S+
	 \t
	 \S+
	 \t
	 (.*)
	 $/x) {
	
	my $start_2 = $1;
	my $end_2 = $2;
	my $strand_2 = $3;
	my $id_2 = $5;
	my $desc_2 = $6;
	$operons->add_vertex($id_2);

	#warn "start: $start_2\n";
	#warn "end_2: $end_2\n";
	#warn "id: $id_2\n";
	#warn "desc: $desc_2\n";
	
	
	$gene_start{$id_2} = $start_2;
	$gene_end{$id_2} = $end_2;
	$id2desc{$id_2} = $desc_2;	
	$id2strand{$id_2} = $strand_2;
	
	if ($strand_1 eq '+' and 
	    $strand_2 eq'+') {
	    
	    ### This region is upstream of gene2
	    $upstream_start{$id_2} = $end_1;
	    $upstream_end{$id_2} = $start_2 + $offset;

	    if ($start_2 - $end_1 < $operon_threshold) {
		$operons->add_edge($id_1, $id_2);
	    }
	}

	
	elsif ($strand_1 eq '-' and 
	    $strand_2 eq'-') {
	    ### This region is upstream of gene1
	    $upstream_start{$id_1} = $end_1 - $offset;
	    $upstream_end{$id_1} = $start_2;
	    if ($start_2 - $end_1 < $operon_threshold) {
		$operons->add_edge($id_2, $id_1);
	    }
	}
	
	elsif ($strand_1 eq '-' and 
	    $strand_2 eq'+') {
	    ### This site is upstream of gene1 and gene2
	    $upstream_start{$id_2} = $end_1 - $offset;
	    $upstream_end{$id_2} = $start_2 + $offset;
	
	    $upstream_start{$id_1} = $end_1 - $offset;
	    $upstream_end{$id_1} = $start_2 + $offset ;
	}
	
	else {
	    ### not upstream of any gene
	}


	($start_1, $end_1, $strand_1, $id_1, $desc_1) = ($start_2, $end_2, $strand_2, $id_2, $desc_2);


    }
}
    


### Calculate frequency distribution of each description word over all operons
my %word2freq;
foreach my $id ($operons->vertices) {
    my $_id = $id;
    my @successors;
    my $first_id = $_id;
    while ( my ($successor) = $operons->successors($_id)) {
	push @successors, $successor;
	$_id = $successor;
    }
    my $desc = "";
    foreach my $_id ($first_id, @successors) {
	$desc .= " $id2desc{$id}";
    }
    my %words;
    foreach my $word ( $desc =~ m/(\w{3,})/g  ) {
	$words{lc$word} = 1
    }
    foreach my $word (keys %words) {
	$word2freq{$word} += 1 / ($operons->vertices);
    }
}



foreach my $hmmsearch_file (@hmmsearch_files) {
    
    my %hmmsearch2gene2hits;
    
    ### Read the hmmsearch results
    
    my $minimum_score = $default_threshold;
    
    warn "parsing file '$hmmsearch_file' using threshold score $minimum_score\n";
	
    my %gene2hits;
    open (FILE, "<$hmmsearch_file") or die "Failed to open '$hmmsearch_file'\n$!";
    while (<FILE>) {
	if (m/
	    ([\d\.]+)
	    \s+
	    \(bits\)\s+
	    f\:
	    \s*
	    (\d+)
	    \s+
	    t\:
	    \s*  
	    (\d+)
	    \s+
	    Target\:
	    /x and
	    $1 > $minimum_score
	    ) {
	    #warn " $1 > $minimum_score\n";
	    
	    my ($score,$start, $end) = ($1,$2,$3);
	    my $strand = "?";
	    if ($start < $end) {
		$strand = '+';
	    }
	    elsif ($start > $end) {
		$strand = "-";
	    } else {
		die;
	    }
	    
	    
	    ### Create a Hit object for this HMMER hit
	    my $hit = new Hit;
	    $hit->score($score);
	    $hit->start($start);
	    $hit->end($end);
	    $hit->strand($strand);
	    
	    ### Read the next few lines, which should contain the nucleotide sequence
	    my $seqline;
	    foreach my $i (1 .. 4) {
		$seqline = <FILE>;
		chomp $seqline;
	    }
	    #warn "$seqline\n";
	    if ($seqline =~ m/(\d+)\s+(\S+)\s+(\d+)/) {
		my $seq = $2;
		#warn "$seq ($strand)\n";
		$hit->sequence($seq);
		
	    } else {
		die "Could not get sequence from '$seqline' ($start,$end)";
	    }
	    
	    ### Find which gene(s) this site is upstream of
	    my @genes;
	    my $in_upstream_region = 0;
	    foreach my $id (keys %id2desc) {
		
		if (defined $upstream_start{$id} and
		    defined $upstream_end{$id}) {
		    
		    my $upstream_start = $upstream_start{$id};
		    my $upstream_end = $upstream_end{$id};
		    
		    
		    ### Allow for incorrect gene start-sites
		    #if ($id2strand{$id} eq '+') {
		    #	$upstream_end += $offset;
		    #$upstream_start -= $offset;
		    
		    #    } elsif ($id2strand{$id} eq '-') {
		    #	$upstream_start -= $offset;
		    #$upstream_end += $offset;
		    #} else {
		    #die "Strand = '$id2strand{$id}'";
		    #}
		    
		    if (
			($hit->start > $upstream_start and
			 $hit->start < $upstream_end) or
			
			($hit->end > $upstream_start and
			 $hit->end < $upstream_end)
			) {
			
			
			### Check relative orientations
			$in_upstream_region = 1;
			push @{$gene2hits{$id}}, $hit;
			
		    }
		    
		}
	    }
	}
    }
    
    $hmmsearch2gene2hits{$hmmsearch_file} = \%gene2hits;
    
    
    
    ### Now sort the genes by descending score
    my %score2genes;
    
    foreach my $id (keys %gene2hits) {
	my $best_score = -99;
	foreach my $hit (@{$gene2hits{$id}}) {
	    $best_score = $hit->score if $hit->score > $best_score;
	}	
	$score2genes{$best_score}{$id} = 1;
	
    }
    my @genes;
    foreach my $score (sort {$b<=>$a} keys %score2genes) {
	push @genes, sort keys %{$score2genes{$score}};
    }
    


 
    ### Which words are associated with which hits?
    my %word2hits;
    foreach my $id (@genes) {
	foreach my $word ( $id2desc{$id} =~ m/(\w{3,})/gi ) {
	    $word = lc $word;
	    foreach my $hit ( @{$gene2hits{$id}}) {
		my $start = $hit->start();
		my $end = $hit->end();
		my $strand = $hit->strand();
		$word2hits{$word}{"$strand:$start-$end"} = 1;
		#warn "\$word2hits{$word}{\"$strand:$start-$end\"} = 1\n";
	    }
	}
    }
   
    
    ### How many different hits are there for this motif?
    my %hits;
    foreach my $id (keys %id2desc) {
	foreach my $hit ( @{$gene2hits{$id}}) {
	    my $start = $hit->start();
	    my $end = $hit->end();
	    my $strand = $hit->strand();
	    $hits{"$strand:$start-$end"} = $1;
	}
    }

    my %word2freq_for_this_hmmsearch;


    ### Are there any significantly over-represented words?
    my %significant_words;
    my $has_significant_words=0;
    foreach my $word (sort keys %word2hits) {
	my $n = keys %hits;
	my $p = $word2freq{$word};
	my $expect = $p * $n;
	my $k = keys %{$word2hits{$word}};
	my $significance = get_P_value($n,$k,$p);
	if ($significance < $significance_threshold and $k > 1) {
	    $significant_words{$word} = $significance;
	    $has_significant_words=1;
	}

	### We don't want to keep motifs that are associated with transposons etc
	foreach my $word (keys %significant_words) {
	    $has_significant_words=0 if $word =~ m/transpos|integrase|ISPsy|isxac3/i;
	}
	



    }


    if ($has_significant_words or 1) {
	print "<h3>$hmmsearch_file</h3>\n";
	print "<ul>\n";
	foreach my $word (sort keys %significant_words) {
	    my $n = keys %hits;
	    my $p = $word2freq{$word};
	    my $expect = $p * $n;
	    my $k = keys %{$word2hits{$word}};
	    my @hits = keys %{$word2hits{$word}};
	    #my $significance = get_P_value($n,$k,$p);
	    my $significance = $significant_words{$word};
	    print "<li>'<b>$word</b>' occurs in $k hits out of $n (expect $expect) (P=$significance) (@hits)</li>\n";
	}
	print "</ul>\n";
	
	### Print information for each gene
	print  "\n";
	print  "<table border=1>\n";
	print  "<tr>";
	print  "<th>Target operon</th>";
	print  "<th>Hit sites</th>";
	print  "</tr>";
	
	foreach my $id (@genes) {
	    print  "<tr>";
	    
	    ### print  gene and its downstream genes in same operon
	    print  "<td>";
	    my @successors;
	    my $first_id = $id;
	    while ( my ($successor) = $operons->successors($id)) {
		push @successors, $successor;
		$id = $successor;
	    }
	    foreach my $id ($first_id, @successors) {
		my $desc = $id2desc{$id};
		foreach my $word (keys %significant_words) {
		    $desc =~ s/($word)/<b>$1<\/b>/gi;

		}
		print  "<p><em>$id</em> ";
		print  "($gene_start{$id}-$gene_end{$id} $id2strand{$id}) ";
		print  "$desc";
		print  "</p>\n";
	    }
	    print  "</td>\n";
	    $id = $first_id;
	    
	    
	    ### print  the hits to hmm
	    print  "<td>";
	    foreach my $hit ( @{$gene2hits{$id}}) {
		
		### Does this hit fall inside a gene?
		my @inside_genes = inside_genes(\%gene_start, \%gene_end, $hit->start, $hit->end); 
		
		print  "<p>";
		print  $hit->start()."-".$hit->end()." ";
		print  "(".$hit->strand().") ";
		print  "score=".$hit->score()." ";
		
		my $sequence = $hit->sequence();
		$sequence =~ s/-//g;
		
		print  "<tt>$sequence</tt> ";
		foreach my $gene (sort @inside_genes) {
		    print  " (inside $gene)";	
		}
		print  "</p>";
	    }
	    print  "</td>\n";
	    
	    
	    
	    
	    
	    print  "</tr>\n";
	}
	print  "</table>\n\n";
	
	
	print "<hr></hr>\n";
    }
}




    
print  $cgi->end_html;







sub get_P_value {
  ### Find probability that number of matches is greater than or equal to $k
  ### ie $p will yield >= $k heads in $n flips
  my ($n, $k, $p) = @_ ;
  my $_P = 0;
  for (my $i = 0 ; $i < $k ; $i++) {
    $_P += binomial($n, $i, $p)
  }
  my $P = 1 - $_P ;
  $P = 100 if $p >=  1 ;
  return $P ;
}

sub binomial {
  my ( $n, $k, $p ) = @_ ;
  return $k == 0 if $p == 0 ;
  return $k != $n if $p == 1 ;
  return choose($n, $k) * $p**$k * (1-$p)**($n-$k) ;
}

sub choose {
  my ($n, $k) = @_ ;
  my ($result, $j) = (1, 1) ;
  
  return 0 if $k > $n || $k < 0 ;
  $k = ($n - $k) if ($n - $k) < $k ;

  while ( $j <= $k ) {
    $result *= $n-- ;
    $result /= $j++ ;
  }
  return $result ;
}





sub inside_genes{
    my $gene_start = shift or die;
    my $gene_end = shift or die;
    my $hit_start = shift or die;
    my $hit_end = shift or die;
    my %inside_genes;
    foreach my $gene (keys %{$gene_start}) {
	if ($hit_start >= $$gene_start{$gene} and $hit_start <= $$gene_end{$gene} ) {
	    $inside_genes{$gene} = 1;
	} elsif ($hit_end >= $$gene_start{$gene} and $hit_end <= $$gene_end{$gene} ) {
	    $inside_genes{$gene} = 1;
	}
    }
    return (keys %inside_genes);
}


########################################################################################################################





### Hit.pm class of objects representing a hit,
###  i.e. a good match to the sought sequence
###  defined by the frequency matrix

package Hit ;

sub new {
  my $self = {} ;
  bless $self, "Hit" ;
  return $self ;
}


sub start {
  my $self  =  shift @_ ;
  if ( my $specified = shift @_ ) {
    $$self{START} = $specified 
  }
  return $$self{START}
}

sub end {
  my $self  =  shift @_ ;
  if ( my $specified = shift @_ ) {
    $$self{END} = $specified 
  }
  return $$self{END}
}
sub strand {
  my $self = shift @_ ;
  if ( my $specified = shift @_ ) {
    $$self{STRAND} = $specified 
  }
  return $$self{STRAND}
}


sub score {
  my $self = shift @_ ;
  if ( my $specified = shift @_ ) {
    $$self{SCORE} = $specified 
  }
  return $$self{SCORE}
}



sub sequence {
  my $self = shift @_ ;
  if ( my $specified = shift @_ ) {
    $$self{SEQ} = $specified 
  }
  return $$self{SEQ}
}

sub locus {
  my $self = shift @_ ;
  if ( my $specified = shift @_ ) {
    $$self{LOCUS} = $specified 
  }
  return $$self{LOCUS}
}

sub intergenic_size {
  my $self = shift @_ ;
  if ( my $specified = shift @_ ) {
    $$self{SIZE} = $specified 
  }
  return $$self{SIZE}
}



