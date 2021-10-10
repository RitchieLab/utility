#!/usr/bin/perl

# parses the tabbed output from MDR and
# compiles information about the models stored
# in the file

use strict;
my $VERSION_DATE = '10/15/10';

my $CV = 0;
my $MODELSIZE = 1;

sub usage{
  print "\n\n\t$0 $VERSION_DATE\n";
  print "\n\tusage:\t\t$0 <tab file> <out file>\n\n";
  print "\n\texample:\t$0 27.mdr.tab.txt 27.mdr.sum\n\n";
}

sub check_args{
  if(@ARGV != 2){
    usage();
    exit(1);
  }
  return @ARGV;
}

my %models;
my $last_cv=1;

my ($tabfile, $outfile) = check_args();

# determine largest model size
# my $max_model = 0;
open(IN, $tabfile) or die "$0:  $tabfile:  $!\n\n";
<IN>;
while(<IN>){
  chomp;
  my @info = split/\t/;
  if($info[0] == 1){
    if($info[1] > $MODELSIZE){
      $MODELSIZE = $info[1];
    }
  }
  else{
    last;
  }
}
close(IN);

open(IN, $tabfile) or die "$0: $tabfile:  $!\n\n";
<IN>;
while(<IN>){
  unless(/\d/){
    next;
  }

  chomp;
  my @info = split/\t/;

  my $cv = $info[$CV];
  # set up key for the model
  my @snpids;
  for(my $m=0; $m<$MODELSIZE; $m++){
    push(@snpids, $info[2+$m*3]);
  }

  my $key;
  foreach my $snpid(@snpids){
    $key .= $snpid . '.';
  }

  chop($key);

  # store information for the models
  $models{$key}->{'totalcv'}++;
  $models{$key}->{'totaltrain'} += $info[2+$MODELSIZE*3];
  if($info[2+$MODELSIZE*3] > $models{$key}->{'besttrain'}){
  	$models{$key}->{'besttrain'}=$info[2+$MODELSIZE*3];
  }
  push(@{$models{$key}->{'trainarray'}}, $info[2+$MODELSIZE*3]);
  $models{$key}->{'totaltest'} += $info[3+$MODELSIZE*3];
  push(@{$models{$key}->{'testarray'}}, $info[3+$MODELSIZE*3]);
}

close(IN);

# calculate averages
foreach my $loci(keys %models){
  $models{$loci}->{'avgtrain'} = $models{$loci}->{'totaltrain'}/$models{$loci}->{'totalcv'};
  $models{$loci}->{'avgtest'} = $models{$loci}->{'totaltest'}/$models{$loci}->{'totalcv'};
}


sub by_cv_test{
  $models{$b}->{'besttrain'} <=> $models{$a}->{'besttrain'}
    ||
  $models{$b}->{'avgtrain'} <=> $models{$a}->{'avgtrain'}
    ||
  $models{$b}->{'avgtest'} <=> $models{$a}->{'avgtest'}
    ||
  $models{$b}->{'totalcv'} <=> $models{$a}->{'totalcv'}
}

open(OUT, ">$outfile") or die "$0:  $outfile:  $!\n\n";


for(my $m=1; $m<=$MODELSIZE; $m++){
  print OUT "Locus$m\t";
}
print OUT "Best Training\tAvg Training\tAvg Testing\tCVC\n";


foreach my $loci(sort by_cv_test keys %models){
  my @snps = split/\./, $loci;
  my $currmodsize = 0;



  my $line;
  for(my $m=0; $m<$MODELSIZE; $m++){
    $line .= "$snps[$m]\t";
  }
  print OUT $line;

  print OUT "$models{$loci}->{'besttrain'}";
  print OUT "\t$models{$loci}->{'avgtrain'}";
  print OUT "\t$models{$loci}->{'avgtest'}";
  print OUT "\t$models{$loci}->{'totalcv'}";

  print OUT "\n";
}

close(OUT);

print "Wrote summary file $outfile\n";
