#!/usr/bin/perl
# This is script for extract annotation site for the peak annotation
# command:perl Whole_gene_annotation_list_maker.pl hg38_refGene

# The gene is not overlap with other genes
# The gene should have unique TSS and TES
# exon will be marked, and intron between exon will be marked

open FILE, $ARGV[0];
my $first=<FILE>;
my $inputfile=$ARGV[0];
open ( my $outfile, ">whole_gene_$inputfile" ) or die;
print $outfile "gene chr strand TSS TES varient exon_start exon_end hit\n";
my %gene_starts;
my %gene_ends;
my %transcript_starts;
my %transcript_ends;
my %record;
my %gene_basic;
# reading the gene with unique TSS and TES, and record all TSS and TES of transcripts
while (defined (my $line = <FILE>) ) {
  	chomp $line;#delete the /n of each line
  	my @parameters=split("\t",$line);
  	my $ID=$parameters[1];
    my $chr=$parameters[2];
    my $strand=$parameters[3];
  	my $start=$parameters[4];
  	my $end=$parameters[5];
  	my $gene=$parameters[12];
  	$record{$gene}{$ID}=$line;
  	$transcript_starts{$chr}{$strand}{$ID}=$start;
  	$transcript_ends{$chr}{$strand}{$ID}=$end;
  	if (defined $gene_starts{$gene}) {
  		unless($gene_starts{$gene} eq $start) {
  			$gene_starts{$gene}=0;
  		}
  	}else{
  		$gene_starts{$gene}=$start;
  	}
  	if (defined $gene_ends{$gene}) {
  		unless($gene_ends{$gene} eq $end) {
  			$gene_ends{$gene}=0;
  		}
  	}else{
  		$gene_ends{$gene}=$end;
  	}
  	$gene_basic{$gene}{'chr'}="$chr";
  	$gene_basic{$gene}{'strand'}="$strand";
}

# get the list of not overlaped unique TSS and TES gene list
my @filter_list;
foreach my $gene (sort keys %gene_starts){
	if ($gene_starts{$gene}>0 and $gene_ends{$gene}>0) {#get unique gene
		# get the not overlaped gene
		my $check=0;
		my $chr=$gene_basic{$gene}{'chr'};
		my $strand=$gene_basic{$gene}{'strand'};
		my $start=$gene_starts{$gene};
		my $end=$gene_ends{$gene};
		foreach my $ID (keys %{$transcript_starts{$chr}{$strand}}){
			if ($transcript_starts{$chr}{$strand}{$ID}>$start and $transcript_starts{$chr}{$strand}{$ID}<$end) {#start in gene
				$check++;
			}
			if ($transcript_ends{$chr}{$strand}{$ID}>$start and $transcript_ends{$chr}{$strand}{$ID}<$end) {#end in gene
				$check++;
			}
			if ($transcript_ends{$chr}{$strand}{$ID}>$end and $transcript_starts{$chr}{$strand}{$ID}<$start) {#end in gene
				$check++;
			}
		}
		if ($check eq 0) {
			push @filter_list,$gene;
		}
	}
}

# now marking exons for the gene that have unique TSS and TES, not overlap with other genes
foreach my $gene (@filter_list){
  my $varient=0;
  my %exons;
  # marking the exons 
  foreach my $ID (keys %{$record{$gene}}){
    $varient++;
    my $line=$record{$gene}{$ID};
    my @parameters=split("\t",$line);
    my $exonStarts=$parameters[9];
    my $exonEnds=$parameters[10];
    my @starts=split(",",$exonStarts);
    my @ends=split(",",$exonEnds);
    my $index=0;
    foreach my $start(@starts){
      $exons{$starts[$index]}{$ends[$index]}=$exons{$starts[$index]}{$ends[$index]}+1;
      $index++;
    }
  }
  foreach my $start (sort {$a<=>$b} keys %exons){
    foreach my $end (sort {$a<=>$b} keys %{$exons{$start}}){
      if ($gene_basic{$gene}{'strand'} eq "+") {
        print $outfile "$gene $gene_basic{$gene}{'chr'} $gene_basic{$gene}{'strand'} $gene_starts{$gene} $gene_ends{$gene} $varient $start $end $exons{$start}{$end}\n";
      }else{
        print $outfile "$gene $gene_basic{$gene}{'chr'} $gene_basic{$gene}{'strand'} $gene_ends{$gene} $gene_starts{$gene} $varient $end $start $exons{$start}{$end}\n";
      }
      
    }
  }
}