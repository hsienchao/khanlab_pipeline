#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long qw(GetOptions);
use File::Basename;
use Cwd 'abs_path';

my $a;
my $b = "/data/khanlab/projects/pipeline_production/khanlab_pipeline/ref/gencode_v32_annotation.genes.txt";
my $i =1;
my $j =1;

my $usage = <<__EOUSAGE__;

Usage:

$0 [options]

Options:

  -a  <string>  The main annotation file
  -b  <string>  The file to be appended (default: $b)
  -i  <string>  The ith in a (default: $i)
  -j  <string>  The jth in b (default: $j)
  
__EOUSAGE__



GetOptions (
  'a=s' => \$a,
  'b=s' => \$b,
  'i=i' => \$i,
  'j=i' => \$j
);

if (!$a) {
    die "Please specifiy the main file!\n$usage";
}

open(SEC_FILE, "$b") or die "cannot open file $b"; 
my $header = <SEC_FILE>;
chomp $header;
my $empty_line = "";
my @headers = split(/\t/, $header);
foreach my $h (@headers) {
	$empty_line = $empty_line."\t.";
}
my %sec = ();
while(<SEC_FILE>) {
	chomp;
	my @fields = split(/\t/);
	$sec{$fields[$j-1]} = $_;
}
close(SEC_FILE);


open(MAIN_FILE, "$a") or die "cannot open file $a"; 
my $main_header = <MAIN_FILE>;
chomp $main_header;
print "$main_header\t$header\n";
while(<MAIN_FILE>) {
	chomp;
	my @fields = split(/\t/);
	my $gene_id = $fields[$i-1];
	if (exists $sec{$gene_id}) {
		print "$_\t$sec{$gene_id}\n";
	} else {
		print "$_\t$empty_line\n";
	}	
}
close(MAIN_FILE);
