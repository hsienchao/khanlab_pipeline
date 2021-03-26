#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long qw(GetOptions);
use File::Basename;
use Cwd 'abs_path';

my $sqanti_file;
my $gtf_file;
my $out_file;

my $usage = <<__EOUSAGE__;

Usage:

$0 [options]

Options:

  -c  <string>  Santi3 classification file
  -g  <string>  Input GTF file
  -o  <string>  Filtered GTF file
  
__EOUSAGE__



GetOptions (
  'c=s' => \$sqanti_file,
  'g=s' => \$gtf_file,
  'o=s' => \$out_file  
);

if (!$sqanti_file || !$gtf_file || !$out_file) {
    die "Please specifiy options!\n$usage";
}

open(SQANTI_FILE, "$sqanti_file") or die "cannot open file $sqanti_file"; 
my %trans = ();
<SQANTI_FILE>;
while(<SQANTI_FILE>) {
	chomp;
	my @fields = split(/\t/);
	$trans{$fields[0]} = '';
}
close(SQANTI_FILE);

open(GTF_FILE, "$gtf_file") or die "cannot open file $gtf_file"; 
open(OUT_FILE, ">$out_file") or die "cannot open file $out_file"; 
while(<GTF_FILE>) {
	chomp;
	my ($trans_id) = $_ =~ /.*\stranscript_id "(.*)"/;
	if (exists $trans{$trans_id}) {
		print OUT_FILE "$_\n";
	}
}
close(GTF_FILE);
close(OUT_FILE);
