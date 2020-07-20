#!/usr/bin/perl
# This is script for extract bed from mRNA list download from UCSC
# command:perl UCSC2bed.pl mRNA

open FILE, $ARGV[0];
my $first=<FILE>;
open ( my $out, ">hg38_mRNA.bed" ) or die;
while (defined (my $line = <FILE>) ) {
  	chomp $line;#delete the /n of each line
  	my @parameters=split("\t",$line);
  	my $ID=$parameters[1];
    my $chr=$parameters[2];
    my $strand=$parameters[3];
    my $start=$parameters[4];
    my $end=$parameters[5];
    my $gene=$parameters[12];
    print $out "$chr\t$start\t$end\t$gene\t$ID\t$strand\n";
}
