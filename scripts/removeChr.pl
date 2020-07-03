#!/usr/bin/perl -w
use strict;
use Getopt::Long qw(GetOptions);
use Time::Piece;
use File::Basename;
use Cwd;

my $input_file;
my $output_file;
my $chr_list_file;

my $usage = <<__EOUSAGE__;

Usage:

$0 [options]

required options:

	-i <string> Input BED file  
	-o <string> Output BED file
	-l <string> Chromosome list file
	
__EOUSAGE__



GetOptions (
	'i=s' => \$input_file,
	'o=s' => \$output_file,
	'l=s' => \$chr_list_file
);

unless ($input_file && $output_file && $chr_list_file) {
	die "Please input, output file and scale factor\n$usage";
}

my %exclude_chr_list = ();
&getExcludeChrList($chr_list_file);
&removeExcludeChrList($input_file,$output_file);

sub getExcludeChrList {
    my ($exlude_chr_file) = @_;
    open(FILE, $exlude_chr_file) or die "Cannot open file $exlude_chr_file";
    while(<FILE>){
         chomp;
		 $exclude_chr_list{$_} = '';
	}
	close(FILE);
}

sub removeExcludeChrList {
    my ($in_file, $out_file) = @_;
    open(IN_FILE, $in_file) or die "Cannot open file $in_file";
	open(OUT_FILE, ">$out_file") or die "Cannot open file $out_file";
    while(<IN_FILE>){    
         chomp;
		 my @fields = split(/\t/);
		 if (!exists $exclude_chr_list{$fields[0]}) {
		 	print OUT_FILE $_."\n";
		 }
		 
	}
	close(IN_FILE);
	close(OUT_FILE);
}