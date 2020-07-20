#!/usr/bin/perl
# This script is use for cut sam file into centain length

open FILE, $ARGV[0];
my $header=0;

while (defined (my $line = <FILE>) ) {
  	chomp $line;#delete the /n of each line
    my @parameters=split(/\t/,$line);
    if ($parameters[0]=~/^@/) {
    	$header++;
    }else{
    	last;
   	}
}
close FILE;


open FILE, $ARGV[0];
my $cut=$ARGV[1]+$header;
print "The input file have $header lines header, cutting the first $cut lines of file\n";
my $filename=$ARGV[0];
open ( my $outfile, ">cut\_$filename" ) or die; #data output
my $index=0;
while (defined (my $line = <FILE>) ) {
	$index++;
	print $outfile "$line";
	if ($index eq $cut) {
		last;
	}
}
