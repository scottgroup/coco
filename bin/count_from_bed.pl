#!/usr/bin/perl

#
# This script is designed to be used in combination with others, see https://github.com/mw55309/RNAfreak
#
# The script takes a SAM output file from htseq-count, and pipes the non-uniquely mapped
# reads through samtools and BEDTools to create "multi-map groups" - that is, groups of
# genes that multi-mapped reads uniquely map to.  It then counts the number of MMGs 
# created and outputs to STDOUT
#
# Run without any arguments to see usage
#
# NB if you have not used Ensembl genes in the pipeline, you need to change "XF:Z:ENS" below
#

use strict;
use POSIX;


# we should have one argument, the output from featureCounts
unless (@ARGV==1) {
	warn "Usage: perl count_from_bed.pl <output_from_featureCounts> \n";
	exit;
}

# get BED file from the command line
my $bed = shift;
unless (-f $bed) {
	warn "$bed is not a file\n";
	exit;
}


open(IN, "grep  'Assigned' $bed|");


# variables to hold the results
my $r = undef;
my $g - undef;

while(<IN>) {

	# split on whitespace
	my @d = split(/\t+/);

	# get the read ID
	my $read = $d[0];

	# get the hit ID
	my $gene = $d[2];
	$gene =~ s/\n//g;
	
	# count that gene for that read
	# only once using hash refs
	$g->{$read}->{$gene}++;
}
close IN;

# go through $g and count
# instances of gene groups
my %c;
while(my($read,$hr) = each %{$g}) {
	
	# get gene list for this read
	my @genes = sort keys %{$g->{$read}};

	# create a key
	my $genekey = join(",", @genes);

	# count that key in a hash
	$c{$genekey}++;
}

# print out group counts
while(my($group,$count) = each %c) {
	print "$group\t$count\n";
}	

