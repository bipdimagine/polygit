#!/usr/bin/env perl
use strict;
use warnings;

while (<STDIN>) {
    chomp;

    # garder le header tel quel
    if (/^#/) {
        print "$_\n";
        next;
    }

    my @f = split(/\t/);

    # FORMAT + sample (mono-sample)
    my @fmt = split(/:/, $f[8]);
    my @val = split(/:/, $f[9]);

    my (%idx, $i);
    for ($i = 0; $i < @fmt; $i++) {
        $idx{$fmt[$i]} = $i;
    }

    # besoin de GT et AD
    unless (exists $idx{GT} && exists $idx{AD}) {
        print "$_\n";
        next;
    }

    my $gt = $val[ $idx{GT} ];
    my $ad = $val[ $idx{AD} ];

    # AD attendu : ref,alt
    my ($ref, $alt) = split(/,/, $ad);
    if (!defined $ref || !defined $alt || ($ref + $alt) == 0) {
        print "$_\n";
        next;
    }

    my $vaf = $alt / ($ref + $alt);

    # règle de correction
    if ($gt eq '1/1' && $vaf < 0.8) {
        $val[ $idx{GT} ] = '0/1';
        $f[9] = join(":", @val);
    }

    print join("\t", @f), "\n";
}