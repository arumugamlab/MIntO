#!/usr/bin/perl

# '''
# compiles annotation for dbCAN based on overview.txt
# uses all annotations, not just ones with 2 or 3 hits among diamond, hotpep and hmmer

# Authors: Vithiagaran Gunalan
# '''
use strict;

open FILE, $ARGV[0];
my $header = <FILE>;

print "ID\tdbCAN.mod\tdbCAN.enzclass\n";
while (<FILE>){
	my $line = $_;
	$line =~ s/\t\-/\t/g;
	my @array = split /\t/, $line;
	my $id =  shift @array;
	pop @array;
	my $enz;
	my (%mod,%enzymes);
	my @filler = ("-","-");
	my $len = scalar(@array);
	for (my $i=1;$i<=$len;$i++){
		my @hitcols = split /\+/, $array[$i];
		foreach my $hit (@hitcols){
			#$hit =~ s/\-//g;
			#print $hit, "\n";
			if ($hit=~ /(.*?)\(/){
				$enz = $1;
			} else {
				$enz = $hit;
			}
			if ($enz =~ /^CBM/){
				#push @mod, $enz;
        	                $mod{$enz} = "";
                #	} elsif ($entry !~ /^CBM/ && $entry =~ /^\D/){
                	} else {
                        	#push @enzymes, $enz;
				$enzymes{$enz} = "";
			}
		}
	}
	if (%mod){
		my @mod = keys %mod;
                $filler[0] = join ",", @mod;
        }
        if (%enzymes){
        	my @enzymes = keys %enzymes;
                $filler[1] = join ",", @enzymes;
        }
	my $myline =  $id."\t".join "\t", @filler;
        print $myline, "\n";
}
