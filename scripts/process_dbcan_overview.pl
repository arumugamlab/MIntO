#!/usr/bin/env perl

# '''
# compiles annotation for dbCAN based on overview.txt
# uses all annotations, not just ones with 2 or 3 hits among diamond, hotpep and hmmer

# Authors: Vithiagaran Gunalan, Mani Arumugam
# '''
use strict;
use warnings;

open FILE, $ARGV[0];
my $header = <FILE>;
chomp($header);
my @fields = split(/\t/, $header);

print "ID\tdbCAN.module\tdbCAN.enzclass\tdbCAN.subfamily\teCAMI.subfamily\teCAMI.submodule\tdbCAN.EC\n";
while (<FILE>){
	my $line = $_;
	$line =~ s/\t\-/\t/g;
	my @array = split /\t/, $line;
	my $len = scalar(@array);
    my ($id, $ec);
	my @filler = ("-")x6;
    my $enz;
    my (%mod,%enzymes, %subfams, %ecamisub, %ecamisubmod);

	for (my $i=0;$i<$len;$i++){
        if ($fields[$i] eq "Gene ID") {
            $id =  $array[$i];
        }
        if ($fields[$i] eq "EC#") {
            $ec = $array[$i];
            if ($ec) {
                my @ec = split('\|', $ec);
                my %ec = map { (split(':', $_))[0] => 1 } @ec;
                $ec = join(",", keys(%ec));
            }
        }
        if ($fields[$i] =~ /HMMER|DIAMOND|dbCAN_sub|HotPep/) {
            my @hitcols = split /\+/, $array[$i];
            foreach my $hit (@hitcols){
                #$hit =~ s/\-//g;
                #print $hit, "\n";
                if ($hit=~ /(.*?)\(/){
                    $enz = $1;
                } else {
                    $enz = $hit;
                }
                my ($fam, $subfam) = split('_', $enz);
                if ($fam =~ /^CBM/){
                    $mod{$fam} = "";
                } else {
                    $enzymes{$fam} = "";
                }
                if ($subfam) {
                    if ($fields[$i] eq "dbCAN_sub") {
                        if ($fam =~ /^CBM/){
                            $ecamisubmod{$enz} = "";
                        } else {
                            $ecamisub{$enz} = "";
                        }
                    } else {
                        $subfams{$enz} = "";
                    }
                }
            }
        }
	}
	if (%mod){
		my @mod = sort {$a cmp $b} keys %mod;
		$filler[0] = join ",", @mod;
	}
	if (%enzymes){
		my @enzymes = sort {$a cmp $b} keys %enzymes;
		$filler[1] = join ",", @enzymes;
	}
	if (%subfams){
		my @subfams = sort {$a cmp $b} keys %subfams;
		$filler[2] = join ",", @subfams;
	}
	if (%ecamisub){
		my @ecamisub = sort {$a cmp $b} keys %ecamisub;
		$filler[3] = join ",", @ecamisub;
	}
	if (%ecamisubmod){
		my @ecamisubmod = sort {$a cmp $b} keys %ecamisubmod;
		$filler[4] = join ",", @ecamisubmod;
	}
    if ($ec) {
        $filler[5] = $ec;
    }
	my $myline =  $id."\t".join "\t", @filler;
	print $myline, "\n";
}
