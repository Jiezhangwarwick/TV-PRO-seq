#!/usr/bin/perl
# This script is use for merging bedgraph files of PRO-seq results

# file will be: chr position total file1..fineN
# example: perl TV_bedGraph_merger.pl index1_m.bedGraph index4_m.bedGraph index56_m.bedGraph index8_m.bedGraph


# checking input files
my $file_number=$#ARGV+1;
print "Have $file_number input files\n";
for (my $var = 0; $var < $file_number; $var++) {
	unless (-e $ARGV[$var]) {
		$current_file=$var+1;
		die "Cannot open file$current_file $ARGV[$var]\n";
	}
}

# loading input data (print title and sencond line as total reads)
my %merge;
my @total_reads;
open ( my $outfile, ">Merged.bedGraph" ) or die; #data output
print $outfile "chr position total_reads";

for (my $var = 0; $var < $file_number; $var++) {
	my $current_file_number=$var+1;
	print $outfile " $ARGV[$var]";
	open FILE, $ARGV[$var];
	print "Reading data from $ARGV[$var]\n";
	my $current_total=0;
	while (defined (my $line = <FILE>) ) {
		chomp $line;#delete the /n of each line
  	    my @parameters=split("\t",$line);
  	    my $chr=$parameters[0];
  	    my $start=$parameters[1];
  	    my $end=$parameters[2];
  	    my $reads=$parameters[3];
  	    unless ($reads eq 0) {
  	    	for (my $i = $start; $i <$end ; $i++) {
  	    		$merge{$chr}{$i+1}{$var}=$reads;
  	    		$current_total=$current_total+$reads;
  	    	}
  	    }
  	    
	}
	push @total_reads,$current_total;
}
print $outfile "\n";
print $outfile "total reads perfile";
my @empty_reads;#creat a empty array with input files number for loading reads per file
for (my $var = 0; $var < $file_number; $var++) {
	print $outfile " $total_reads[$var]";
	push @empty_reads,0;
}
print $outfile "\n";
print "Sorting input reads\n";
foreach my $chr (sort keys %merge){
	my $loc_sum=scalar keys %{$merge{$chr}};
	if($loc_sum>10000){
		print "Now sorting reads of $chr\n";
	}
	my %pos_sort;
	foreach my $pos (sort keys %{$merge{$chr}}){
		my $sum_reads=0;
		my @reads_record=@empty_reads;
		foreach my $input_number( sort{{$a} <=> {$b}} keys %{$merge{$chr}{$pos}}){
			my $read=$merge{$chr}{$pos}{$input_number};
			@reads_record[$input_number]=$read;
			$sum_reads=$sum_reads+$read;
		}
		$pos_sort{$pos}= "$chr $pos $sum_reads @reads_record\n";
	}
	foreach my $pos (sort {$a<=>$b} keys %pos_sort){
		print $outfile "$pos_sort{$pos}";
	}
}





