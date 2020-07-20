#!/usr/bin/perl
# This script is use for annotated histone CHIP-seq data to TSS
# bedgraph must be bga format
# perl Liner_bedgraph_v4.pl H3K27ac_input_2.bedgraph UniTSS_mRNA
# perl Liner_bedgraph_v4.pl H3K27ac_input_2.bedgraph Beta_summary

# Have three part of script:
# 1. loading annotation list into hash
# 2. Reading bedgraph and recording the location of reads location and recording peaks when meet the peak location
# 2plus. finish the output and move to next peak, if the region have overlap with last one, record part of it.

############################### part1 ##################################
# 1.1 detect the input file type
open FILE, $ARGV[1];
print "Start loading annotation\n";
$filename=$ARGV[0];

my $title=<FILE>;
chomp $title;
my @titles=split(" ",$title);
my $number=0;
my $chr_title_number;
foreach my $title_name (@titles){
	if ($title_name eq "strand") {
		$strand_title_number=$number;
	}
	$number++;
}

# now, peak file chr loc is 0 and TSS file is 1
if ($strand_title_number eq 1) {
	$filename="TVPeak_$filename";
}else{
  $filename="TSS_$filename";
}
open ( my $outfile, ">$filename" ) or die; #data output
print $outfile "@titles";
for (my $var = -1000; $var < 1001; $var++) {
  print $outfile " $var";
}
print $outfile "\n";


# 1.2 loading annotations:
my %TSS_list;
while (defined (my $line = <FILE>) ) {
  chomp $line;#delete the /n of each line
  if($strand_title_number eq 1){
    my @parameters=split(" ",$line);
    my $chr=$parameters[0];
    my $strand=$parameters[1];
    my $TSS=$parameters[2];
    $TSS_list{$chr}{$TSS}{$strand}=$line;
  }else{
    my @parameters=split("\t",$line);
    my $chr=$parameters[0];
    my $strand=$parameters[2];
    my $TSS=$parameters[1];
    $TSS_list{$chr}{$TSS}{$strand}="@parameters";

  }
}



############################### part2 ##################################
# 2.1 reading the bedgraph file
open FILE, $ARGV[0];
print "Start loading bedgraph\n";
my $last_chr="none";
my %reads_index;

my @empty_array;
for (my $var = 0; $var < 2002; $var++) {
  push @empty_array,0;
}
my @starts=@empty_array;
my @ends=@empty_array;
my @reads=@empty_array;
my @TSS_loc;
my @TSS_information;
my @TSS_strand;
my $scan_number=0;
my $chr_gene_number=0;

while (defined (my $line = <FILE>) ) {
    chomp $line;#delete the /n of each line
    my @parameters=split("\t",$line);
    if($parameters[0] eq $last_chr){
      push @starts,$parameters[1];
      push @ends,$parameters[2];
      push @reads,$parameters[3];
      shift @starts;
      shift @ends;
      shift @reads;
      if ($scan_number eq $chr_gene_number) {
        next;
      }
      # when the new input line end be later than TSS+999, it means all the TSS +-1000 region in the window:
      if ($parameters[2]>$TSS_loc[$scan_number]+999) {
        for (my $var = 0; $var < 2002; $var++) {
            # first line:
            if ($starts[$var]<$TSS_loc[$scan_number]-1000 and $ends[$var]>$TSS_loc[$scan_number]-1001) {
              
              if ($ends[$var]>$TSS_loc[$scan_number]+999) {
                @out_array=($reads[$var])x(2001);
                my $out_number=@out_array;
              }else{
                my $start_number=$ends[$var]-($TSS_loc[$scan_number]-1001);
                @out_array=($reads[$var])x($start_number);
              }
            }
            # middle line:
            if ($starts[$var]>$TSS_loc[$scan_number]-1001 and $ends[$var]<$TSS_loc[$scan_number]+1000) {
              @out_array=(@out_array,($reads[$var])x($ends[$var]-$starts[$var])); 
            }
            # last line;
            if ($ends[$var]>$TSS_loc[$scan_number]+999 and $starts[$var]>$TSS_loc[$scan_number]-1000) {
              my $end_number=$TSS_loc[$scan_number]+1000-$starts[$var];
              @out_array=(@out_array,($reads[$var])x($end_number));
              my $out_number=@out_array;
              my $test=$ends[$var]-$starts[$var];
              last;
            }

        }
        #check if the @out_array is right:
          my $out_number=@out_array;

          if ($out_number eq 2001 and $parameters[0] eq $last_chr) {
            if ($TSS_strand[$scan_number] eq '+') {
              print $outfile "$TSS_information[$scan_number] @out_array\n";
            }else{
              @out_array=reverse @out_array;
              print $outfile "$TSS_information[$scan_number] @out_array\n";
            }
          }else{
            if($out_number eq 0 and $parameters[0] eq $last_chr){
              @out_array=('0')x(2001);
              print $outfile "$TSS_information[$scan_number] @out_array\n";
            }else{
              #print "Following line's number is $out_number\nAnnotation is $TSS_information[$scan_number]\nlast line is $line\n";
            }  
          }

          
          $scan_number++;
          if ($scan_number%1000 eq 0){
            print "print $scan_number target for $last_chr\n";
          }

      }
    }else{
      # print out the gene have not been print out:
      if ($scan_number<$chr_gene_number) {
        for (my $gene_index = $scan_number; $gene_index < $chr_gene_number; $gene_index++) {
          my @out_array;
          for (my $var = 0; $var < 2002; $var++) {
            # first line:
            if ($starts[$var]<$TSS_loc[$gene_index]-1000 and $ends[$var]>$TSS_loc[$gene_index]-1000) {

              if ($ends[$var]>$TSS_loc[$gene_index]+999) {
                @out_array=($reads[$var])x(2001);
              }else{
                my $start_number=$ends[$var]-($TSS_loc[$gene_index]-1001);
                @out_array=($reads[$var])x($start_number);
              }
            }
            # middle line:
            if ($starts[$var]>$TSS_loc[$gene_index]-1001 and $ends[$var]<$TSS_loc[$gene_index]+1000) {
              @out_array=(@out_array,($reads[$var])x($ends[$var]-$starts[$var]));
            }
            # last line;
            if ($ends[$var]>$TSS_loc[$gene_index]+999) {
              @out_array=(@out_array,($reads[$var])x(($TSS_loc[$gene_index]+1000)-$starts[$var]));
              last;
            }
          }
          #check if the @out_array is right:
          my $out_number=@out_array;
          if ($out_number eq 2001 and $parameters[0] eq $last_chr) {
            if ($TSS_strand[$scan_number] eq '+') {
              print $outfile "$TSS_information[$scan_number] @out_array\n";
            }else{
              @out_array=reverse @out_array;
              print $outfile "$TSS_information[$scan_number] @out_array\n";
            }
          }else{
            if($out_number eq 0 and $parameters[0] eq $last_chr){
              @out_array=('0')x(2001);
              print $outfile "$TSS_information[$scan_number] @out_array\n";
            }else{
              #print "Following line's number is $out_number\nAnnotation is $TSS_information[$scan_number]\nlast line is $line\n";
            }   
          }
          
        }

      }

      # loading in the TSS locations and reset the arrays
      $chr_gene_number=0;
      $last_chr=$parameters[0];
      undef @TSS_loc;
      undef @TSS_information;
      undef @TSS_strand;
      foreach my $TSS (sort {$a<=>$b} keys %{$TSS_list{$last_chr}}){
        foreach my $strand (sort keys %{$TSS_list{$last_chr}{$TSS}}){
          push @TSS_loc,$TSS;
          push @TSS_information,$TSS_list{$last_chr}{$TSS}{$strand};
          push @TSS_strand,$strand;
          $chr_gene_number++;
        }
      }
      if($chr_gene_number>100){
        print "start counting $last_chr, have $chr_gene_number target\n";
      }
      $scan_number=0;
      @starts=@empty_array;
      @ends=@empty_array;
      @empty_array=@reads;
      
      push @starts,$parameters[1];
      push @ends,$parameters[2];
      push @reads,$parameters[3];
      shift @starts;
      shift @ends;
      shift @reads;
    }
}










