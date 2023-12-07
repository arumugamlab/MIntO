#!/usr/bin/env perl
use strict;
use warnings;

my $Desc = {};
while (<>) {
    if (m/^Funct/) {
        print $_;
        next;
    }
    chomp();
    my ($entry, $desc) = split(/\t/, $_);
    $Desc->{$entry} = $desc;
}
foreach my $entry (sort {$a cmp $b} keys %$Desc) {
    printf "%s\t%s\n", $entry, $Desc->{$entry};
}
