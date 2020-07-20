#!/usr/bin/perl
# This is script annotated peaks with in annotation list:

# command:perl rFAM_region.pl rFAM_annotation Beta_summary

open FILE, $ARGV[0];
my %annotation;
print "loading annotation table\n";

while (defined (my $line = <FILE>) ) {
  	chomp $line;#delete the /n of each line
  	my @parameters=split("\t",$line);
    my $chr=$parameters[0];
    my $strand=$parameters[1];
    my $start=$parameters[2];
    my $end=$parameters[3];

    $annotation{$chr}{$strand}{$start}{$end}="$parameters[4]\t$parameters[5]";
} 
close FILE;


open FILE, $ARGV[1];
my $first=<FILE>;
my $peaks;
print "loading peaks table\n";

while (defined (my $line = <FILE>) ) {
  	chomp $line;#delete the /n of each line
  	my @parameters=split(" ",$line);
    my $chr=$parameters[0];
    my $strand=$parameters[1];
    my $pos=$parameters[2];
    my $beta=$parameters[7];
    $peaks{$chr}{$strand}{$pos}=$beta;
} 
close FILE;

my $filename=$ARGV[0];


open ( my $outfile, ">rFAM_beta" ) or die;
print $outfile "chr\tstrand\tpeak\tbeta\tclass\ttype\tstart\tend\n";
print "annotating peaks table\n";
foreach my $chr (sort keys %peaks){
	foreach my $strand (sort keys %{$peaks{$chr}}){
		#loading the annotation for the strand:
		my @start_list;
		my @end_list;
		my @name_list;
		foreach my $start (sort {$a<=>$b} %{$annotation{$chr}{$strand}}){
			foreach my $end (sort {$a<=>$b} %{$annotation{$chr}{$strand}{$start}}){
				push @start_list,$start;
				push @end_list,$end;
				push @name_list,$annotation{$chr}{$strand}{$start}{$end};
			}
		}
		my $length=@name_list;
		if ($length>5000) {
			print "Now annotation for $chr\n";
		}
		my $index=0;
		foreach my $pos (sort {$a<=>$b} keys %{$peaks{$chr}{$strand}}){
			while ($pos>$end_list[$index]) {
				if ($index eq $length) {
					last;
				}
				$index++;
			}
			if ($index eq $length) {
				last;
			}
			if ($pos > $start_list[$index]) {
				print $outfile "$chr\t$strand\t$pos\t$peaks{$chr}{$strand}{$pos}\t$name_list[$index]\t$start_list[$index]\t$end_list[$index]\n";
			}
		}
	}
}
close $outfile;

