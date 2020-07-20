#!/usr/bin/perl
# This script is use for change the information of sam file to enlong the location to the end.

# Two things to do:
# 1. change all S to M
# 2. If having S before M, minus the location for the number of S

# Additional function, cut sam file to certain reads:
# The cut number should be put after input file name

# reading the chr length in header:
open FILE, $ARGV[0];
my $header=0;
my %chr_length;
while (defined (my $line = <FILE>) ) {
  	chomp $line;#delete the /n of each line
  	$line=~s/SN://;
  	$line=~s/LN://;
    my @parameters=split(/\t/,$line);
    if ($parameters[0]=~/^@/) {
    	$header++;
    	if ($parameters[0] eq '@SQ') {
    		$chr_length{$parameters[1]}=$parameters[2];

    	}
    }else{
    	last;
   	}
}
print "The hearer have $header lines\n";
close FILE;

open FILE, $ARGV[0];
my $filename=$ARGV[0];
open ( my $outfile, ">Enlong\_$filename" ) or die; #data output
open ( my $report, ">report" ) or die; #data output
print $report "chr\tchr_length\tread_start\tread_end\tmiss_type\tline_number\n";
# directly copy the header:
for (my $var = 0; $var < $header; $var++) {
	my $line = <FILE>;
	print $outfile "$line";
}

# elong the soft fitting end to theory genome location
my $out_counter=0;
my $line_conter=0;
while (defined (my $line = <FILE>) ) {
  	chomp $line;#delete the /n of each line
  	$line_conter++;
    my @parameters=split(/\t/,$line);
    if ($parameters[5]=~/^([0-9]+)S([0-9]+)M/) {
    	my $length=$1+$2;
        my $out="$length";
        $out.="M";
		$parameters[5]=~s/^([0-9]+)S([0-9]+)M/$out/;
		$parameters[3]=$parameters[3]-$1;
	}
    if ($parameters[5]=~/([0-9]+)M([0-9]+)S$/) {
		my $length=$1+$2;
        my $out="$length";
        $out.="M";
		$parameters[5]=~s/([0-9]+)M([0-9]+)S$/$out/;
	}
	if ($parameters[3]>0 and $parameters[3]+$parameters[5]<$chr_length{$parameters[2]}){
		$out_counter++;
		my $outline=join("\t",@parameters);
    	print $outfile "$outline\n";
	}else{
		$end=$parameters[3]+$parameters[5];
		if ($parameters[3]>0) {
			$miss_type="end_overflow";
		}else{
			$miss_type="start_overflow";
		}
		print $report "$parameters[2]\t$chr_length{$parameters[2]}\t$parameters[3]\t$end\t$miss_type\t$line_conter\n";
	}
}
close FILE;
print "Extend of soft fitting reads have been finished, $out_counter reads hace been recorded\n";






