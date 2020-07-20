#!/usr/bin/perl
# This is script for extract annotation site for the peak annotation
# command:perl Unique_annotation_maker.pl hg38_refGene

# This script is use for extract unique TSS/TES/splicing site and alternative TSS/TES/splicing

# Outcome file:
# Outcome column:
# chr loc strand gene type number_min number_max hit variant

# Unique TSS:
# hit=variant, number_min=number_max=1, type=start

# Unique TES:
# hit=variant, number_max=-1, type=end

# Unique splicing site:
# hit=variant, number_min>1, type=end/start

open FILE, $ARGV[0];
my $first=<FILE>;
my $inputfile=$ARGV[0];

my %gene_variant;	# $gene_variant{gene}
my %exon_start;		# $exon_start{gene}{loc}
my %exon_end;		# $exon_end{gene}{loc}
my %gene_basic;		# $gene_basic{gene}{'chr'}=chr	$gene_basic{gene}{'strand'}=strand


while (defined (my $line = <FILE>) ) {
  	chomp $line;#delete the /n of each line
  	my @parameters=split("\t",$line);
  	my $ID=$parameters[1];
    my $chr=$parameters[2];
    my $strand=$parameters[3];
  	my $exonCount=$parameters[8];
  	my $exonStarts=$parameters[9];
  	my $exonEnds=$parameters[10];
  	my $gene=$parameters[12];

  	my @starts=split(",",$exonStarts);
  	my @ends=split(",",$exonEnds);
    if($strand eq "-"){
      my @temp_ends=reverse @starts;
      @starts=reverse @ends;
      @ends=@temp_ends;
    }

    if (defined $gene_variant{$gene}) {
    	$gene_variant{$gene}=$gene_variant{$gene}+1;
    }else{
    	$gene_variant{$gene}=1;
    }

    $gene_basic{$gene}{'chr'}=$chr;
    $gene_basic{$gene}{'strand'}=$strand;

    for (my $i = 0; $i < $exonCount; $i++) {
    	my $number=$i+1;
    	if (defined $exon_start{$gene}{$starts[$i]}) {
        $exon_start{$gene}{$starts[$i]}="$exon_start{$gene}{$starts[$i]} $number";
      }else{
        $exon_start{$gene}{$starts[$i]}=$number;
      }
    	if ($number eq $exonCount) {
    		$number=(-1);
    	}
      if (defined $exon_end{$gene}{$ends[$i]}) {
        $exon_end{$gene}{$ends[$i]}="$exon_end{$gene}{$ends[$i]} $number";
      }else{
        $exon_end{$gene}{$ends[$i]}=$number;
      }
    }
}
close FILE;


open ( my $out, ">All_$inputfile" ) or die;
open ( my $out1, ">UniTSS_$inputfile" ) or die;
open ( my $out2, ">UniTES_$inputfile" ) or die;
open ( my $out3, ">UnSplicing_$inputfile" ) or die;

# chr loc strand gene type number_min number_max hit variant
print $out "chr\tloc\tstrand\tgene\ttype\tnumber_min\tnumber_max\thit\tvariant\n";
print $out1 "chr\tloc\tstrand\tgene\ttype\tnumber_min\tnumber_max\thit\tvariant\n";
print $out2 "chr\tloc\tstrand\tgene\ttype\tnumber_min\tnumber_max\thit\tvariant\n";
print $out3 "chr\tloc\tstrand\tgene\ttype\tnumber_min\tnumber_max\thit\tvariant\n";
# First write the All_$inputfile
foreach my $gene (sort keys %exon_start){
  foreach my $loc (sort {$a<=>$b} keys %{$exon_start{$gene}}){
    my $chr=$gene_basic{$gene}{'chr'};
    my $strand=$gene_basic{$gene}{'strand'};
    my @numbers=split(" ",$exon_start{$gene}{$loc});
    @numbers=sort{$a<=>$b}@numbers;
    $number_min=$numbers[0];
    $number_max=$numbers[-1];
    $hit=@numbers;
    $variant=$gene_variant{$gene};
    print $out "$chr\t$loc\t$strand\t$gene\tstart\t$number_min\t$number_max\t$hit\t$variant\n";
  }
}


foreach my $gene (sort keys %exon_end){
  foreach my $loc (sort {$a<=>$b} keys %{$exon_end{$gene}}){
    my $chr=$gene_basic{$gene}{'chr'};
    my $strand=$gene_basic{$gene}{'strand'};
    my @numbers=split(" ",$exon_end{$gene}{$loc});
    @numbers=sort{$a<=>$b}@numbers;
    $number_min=$numbers[0];
    $number_max=$numbers[-1];
    $hit=@numbers;
    $variant=$gene_variant{$gene};
    print $out "$chr\t$loc\t$strand\t$gene\tend\t$number_min\t$number_max\t$hit\t$variant\n";
  }
}

open FILE,"All_$inputfile";

while (defined (my $line = <FILE>) ) {
    chomp $line;#delete the /n of each line
    my @parameters=split("\t",$line);
    my $type=$parameters[4];
    my $number_min=$parameters[5];
    my $number_max=$parameters[6];
    my $hit=$parameters[7];
    my $variant=$parameters[8];

    if ($number_max eq 1 and $hit eq $variant and $type eq "start") {
      print $out1 "$line\n";
    }
    if ($number_max eq -1 and $hit eq $variant and $type eq "end") {
      print $out2 "$line\n";
    }
    if ($number_min>1 and $hit eq $variant) {
      print $out3 "$line\n";
    }
}




