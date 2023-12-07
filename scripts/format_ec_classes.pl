#! /usr/bin/env perl
use strict;
use warnings;

my $Def = {};

# Get info
while (<>) {

    chomp();

    next unless m/\.- /; # Only definitions

    # Format line
    s/\.\s([\d-])/.$1/g;
    s/\.$//;
    s/\s+/\t/;

    my ($ec, $def) = split(/\t/, $_);
    $Def->{$ec} = $def;
}

# Augment info
my @ecs = sort {$a cmp $b} keys %$Def;
foreach my $ec (@ecs) {
    my @fields = split('\.', $ec);
    for (my $i=0; $i<$#fields; $i++) {
        if ($fields[$i+1] eq "-") {
            $fields[$i] = "-";
        }
    }
    my $parent = join('.', @fields);
    #print "$ec ==> $parent\n";
    if ($Def->{$parent}) {
        $Def->{$ec} = $Def->{$parent}."; ".$Def->{$ec};
    }
}

# Print info
foreach my $ec (@ecs) {
    printf "%s\t%s\n", $ec, $Def->{$ec};
}
