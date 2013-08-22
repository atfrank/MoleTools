#!/usr/bin/env perl

use strict;
use warnings;

open (META,">metadatafile");

for (my $i=1; $i<=72; $i++){
    open(INP,"butane.$i") || die "Cannot open file \"butane.$i\"\n";
    open (GROSS,">grossfield.$i");
    while (<INP>){
	my @t=split/\s+/,$_;
	print GROSS "$t[0] $t[1]\n";
    }
    my $eq=180-($i-1)*5;
    $eq=$eq+360 if ($eq<0);
    print META "grossfield.$i $eq 0.04\n";
    close(GROSS);
    close(INP);
}

close(META);

print "wham P 0.000 359.999 50 0.00001 300 0 metadatafile freefile\n";
