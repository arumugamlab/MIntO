#!/usr/bin/env perl

# '''
# compile annotations from eggnog-mapper

# Authors: Vithiagaran Gunalan, Mani Arumugam
# '''

use strict;
use warnings;

open FILE, $ARGV[0];

my $header = <FILE>; #retrieve and keep header line
print $header;

while(<FILE>){
	my $line = $_;
	$line =~ s/\R//g;
	if ($line =~ /\t$/){
		$line =~ s/\t$//g;
	}
	my @array = split /\t/, $line;

	# eggNOG OGs
	my @OGs = split /\,/, $array[1];
	my @fixed = ();
	foreach my $item(@OGs){
		if ($item =~ /(.*?)\@\d\|root$/){
			push @fixed, $1;
		}
	}
	if (@fixed){
		$array[1] = join ",", @fixed;
	} else {
		$array[1] = "-"
	}

	# KEGG Pathway
	# duplicated entries as ko12345 and map12345 are unified into just map12345
	my $H = {};
	map { s/^ko/map/; $H->{$_} = 1 } split /\,/, $array[3];
	@fixed = sort keys(%$H);
	if (@fixed){
		$array[3] = join ",", @fixed;
	} else {
		$array[3] = "-"
	}

	# Print them all!
	print join "\t", @array, "\n";
}
