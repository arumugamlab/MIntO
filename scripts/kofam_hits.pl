#!/usr/bin/perl

# '''
# # compiles annotation for KofamScan

# Authors: Vithiagaran Gunalan
# '''

use strict;

open FILE, $ARGV[0];

my %hash;
while (<FILE>){
	my $line = $_;
	$line =~ s/\R//g;
	my @array = split "\t", $line;
	if ($array[1] =~ /^[A-Z]/){
		if (exists($hash{$array[0]})){
			$hash{$array[0]} = $hash{$array[0]}.",".$array[1];
		} else {
			$hash{$array[0]} = $array[1];
		}
		#print $array[0], "\t", $array[1], "\n";
		#sleep 1;
		#$hash{$array[0]} = $array[1];
	}
}
print "ID\tkofam_KO\n";
#my $count = 0;
for my $keys (keys %hash){
	print $keys, "\t", $hash{$keys}, "\n";
}
