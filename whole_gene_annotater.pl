# This is script for extract annotation site for the peak annotation
# command:perl whole_gene_annotater.pl whole_gene_hg38_refGene Above_pois_k_value

# The inputs file have all the exon from gene with unique TSS & TES and does not overlap with other genes.

open FILE, $ARGV[0];
my $first_1=<FILE>;

my %genes;
my %exons;
my %basic_gene;
while (defined (my $line = <FILE>) ) {
  	chomp $line;#delete the /n of each line
  	my @parameters=split(" ",$line);
  	my $gene=$parameters[0];
  	my $chr=$parameters[1];
  	my $strand=$parameters[2];
  	my $TSS=$parameters[3];
  	my $TES=$parameters[4];
  	my $exon_start=$parameters[6];
  	my $exon_end=$parameters[7];
  	if ($strand eq "+") {
  		$genes{$chr}{$strand}{$gene}{'start'}=$TSS;
  		$genes{$chr}{$strand}{$gene}{'end'}=$TES;
  		$exons{$gene}{$exon_start}{$exon_end}=$line;
  	}else{
  		$genes{$chr}{$strand}{$gene}{'start'}=$TES;
  		$genes{$chr}{$strand}{$gene}{'end'}=$TSS;
  		$exons{$gene}{$exon_end}{$exon_start}=$line;
  	}
  	$basic_gene{$gene}="$gene $chr $strand $TSS $TES $parameters[5]";#record information for intron 
}

open FILE, $ARGV[1];
my $first_2=<FILE>;
my $filename=$ARGV[0];
open ( my $outfile, ">$filename\_peaks" ) or die;
print $outfile "gene chr strand TSS TES varient region_start region_end hit $first_2";
while (defined (my $line = <FILE>) ) {
	chomp $line;#delete the /n of each line
	my @parameters=split(" ",$line);
  	my $chr=$parameters[0];
  	my $strand=$parameters[1];
  	my $pos=$parameters[2];
  	my $hit_gene="N";
  	foreach my $gene (keys %{$genes{$chr}{$strand}}){
  		if ($genes{$chr}{$strand}{$gene}{'start'}<$pos and $genes{$chr}{$strand}{$gene}{'end'}>$pos) {
  			$hit_gene=$gene;
  			last;
  		}
  	}
  	# looking into the location of peak located in gene
  	unless ($hit_gene eq "N") {
  		my $anno;
  		my $check=0;
  		foreach my $start (keys %{$exons{$hit_gene}}){
  			foreach my $end (keys %{$exons{$hit_gene}{$start}}){
  				if ($start-1<$pos and $end>$pos) {
  					$anno=$exons{$hit_gene}{$start}{$end};
  					$check++;
  				}
  			}
  		}
  		# now the exon have been marked as check>0
  		if($check > 0){
  			print $outfile "$anno $line\n";
  		}else{
  			my $intron_start=0;
	  		my $intron_end=10000000000000;
  			foreach my $start (keys %{$exons{$hit_gene}}){
  				foreach my $end (keys %{$exons{$hit_gene}{$start}}){
  					if ($start<$pos and $start>$intron_start) {
  						$intron_start=$start;
 	 				}
  					if ($end>$pos and $end<$intron_end) {
  						$intron_end=$end;
  					}
	  			}
 		 	}
  			if ($strand eq "+") {
  				print $outfile "$basic_gene{$hit_gene} $intron_start $intron_end 0 $line\n";
	  		}else{
  				print $outfile "$basic_gene{$hit_gene} $intron_end $intron_start 0 $line\n";
  			}
  		}
  	}
}