#!/usr/bin/env perl

# '''
# # compiles annotation for KofamScan

# Authors: Vithiagaran Gunalan, Mani Arumugam
# '''

use strict;
use warnings;
use Getopt::Long;

my ($module_map, $pathway_map);
GetOptions("module-map=s" => \$module_map,
           "pathway-map=s" => \$pathway_map);

if (!$module_map || !$pathway_map) {
    die "Usage: $0 --module=<file> --pathway=<file> <kofam-output-file>";
}

my $KO2Module = {};
open(F, "<$module_map") || die "Cannot open file $module_map: $!";
while (<F>) {
    chomp();
    my ($module, $ko) = split(/\t/);
    $KO2Module->{$ko}->{$module} = 1;
}
close(F);

my $KO2Pathway = {};
open(F, "<$pathway_map") || die "Cannot open file $pathway_map: $!";
while (<F>) {
    chomp();
    my ($pathway, $ko) = split(/\t/);
    $KO2Pathway->{$ko}->{$pathway} = 1;
}
close(F);

my $Hash = {};
while (<>) {
    chomp();
    my ($gene, @annot) = split("\t", $_);
    if (@annot){
        $Hash->{$gene} = join(",", sort {$a cmp $b} @annot);
    }
}
printf "%s\t%s\t%s\t%s\n", "ID", "kofam_KO", "kofam_Module", "kofam_Pathway";
for my $gene (sort {$a cmp $b} keys %$Hash) {
    my $ko_string = $Hash->{$gene};
    my @kos = split(",", $ko_string);

    # Get modules for the KOs
    my $module_string = "-";
    my %Modules  = map {$_ => 1} map {keys(%{$KO2Module->{$_}})} @kos;
    if (%Modules) {
        my @modules  = sort {$a cmp $b} keys(%Modules);
        $module_string = join(",", @modules);
    }

    # Get pathways for the KOs
    my $pathway_string = "-";
    my %Pathways  = map {$_ => 1} map {keys(%{$KO2Pathway->{$_}})} @kos;
    if (%Pathways) {
        my @pathways  = sort {$a cmp $b} keys(%Pathways);
        $pathway_string = join(",", @pathways);
    }
    printf "%s\t%s\t%s\t%s\n", $gene, $ko_string, $module_string, $pathway_string;
}
