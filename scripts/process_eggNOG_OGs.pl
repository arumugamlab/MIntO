#!/usr/bin/perl

# '''
# compile annotations from eggnog-mapper

# Authors: Vithiagaran Gunalan
# '''

use strict;

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
	my @OGs = split /\,/, $array[1];
	my @fixed = ("-");
	foreach my $item(@OGs){
		if ($item =~ /(.*?)\@\d\|root$/){
			my $og = $1;
			push @fixed, $og;
		}
	}
	if (scalar(@fixed) > 1){
		shift @fixed;
		my $joined = join ",", @fixed;
	        splice @array, 1, 1, $joined;
	} else {
		my $joined = join "", @fixed;
		splice @array, 1, 1, $joined;
	}
	print join "\t", @array, "\n";
}
