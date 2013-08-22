#!/usr/bin/env perl

use strict;
use warnings;

my $dir="/green/s3x/slaw/for/afra/butane";
my $k=0.02;

for (my $i=1; $i<=73; $i++){
    open (OUT,">butane.$i");
    my $cmd="awk '{print \$4}' $dir/analyze.$i";
    my @d=`$cmd`;
    my $l=1;
    foreach my $angle (@d){
	chomp $angle;
	$angle=$angle+360 if ($angle < 0);
	$angle=$angle-360 if ($angle > 360);
	my $out="$l\t$angle";
	for (my $j=0; $j<73; $j++){
	    my $e=(180-$j*5);
	    $e=$e+360 if ($e<0);
	    $e=$e-360 if ($e>360);
	    my $dx=$angle-$e;
	    $dx=abs($dx);
	    if ($dx>180){
		$dx=$dx-360;
	    }
	    my $EU=$k*($dx*$dx);
	    $out.="\t$EU";
	}
	print OUT "$out\n";
	$l++;
    }
    close (OUT);
}
