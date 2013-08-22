#!/usr/bin/env perl

use warnings;
use strict;

my @inp;
my $kB=(1.380650424E-23)*(1.44E20); #Boltzmann constant,Conversion from J/K to kcal/mol*K
my $T=300;
my $kBT=$kB*$T;
my $tol=0.00001;
my $bins=10;

while ($#ARGV>=0){
    if ($ARGV[0] eq "-help"){

    }
    elsif ($ARGV[0] eq "-option"){
	#DO SOMETHING
    }
    else{
	push (@inp, shift @ARGV);
    }
}

my $Nwind=$#inp;
my $Niter=100000;
my @nt=();
my @ebf=();
my @ebf2=();
my @fact=();
my @xi=();
my $min=1E10;
my $max=-1E10;
my @hist=();
my @rep=();
my @ebw=();
my $i;
my $j;
my $k;
my $l;
if ($Nwind==-1){
    die "Please provide at least one input file\n";
}

#Read data into memory
#Get $nt[$k]
for ($i=0; $i<=$#inp; $i++){
    open (INP, "$inp[$i]") || die "Cannot open file \"$inp[$i]\"\n";
    $l=0;
    while (<INP>){
	my @s=split/\s+/,$_;
	$xi[$i][$l]=$s[1];
	$min=$xi[$i][$l] if ($xi[$i][$l] < $min);
	$max=$xi[$i][$l] if ($xi[$i][$l] > $max);
	for (my $m=2; $m<=$#s; $m++){
	    $s[$m]=exp(-$s[$m]/$kBT); #Convert potential to exp(-w/kT)
	}
	$ebw[$i][$l]=[splice(@s,2)]; #Splice out energies
	$l++;
    }
    $nt[$i]=$l; #This is one more than the array size
    close(INP);
}

#Get binsize, set up bins for each window
my $binsize=($max-$min)/$bins;
for ($i=0; $i<=$#inp; $i++){
    for (my $b=0; $b<$bins; $b++){
	$hist[$i][$b]=0;
    }
}

#Generate histogram, find representative closest to center of bin
for ($i=0; $i<=$#inp; $i++){
    for ($l=0; $l<$nt[$i]; $l++){
	my $j=1;
	while ($xi[$i][$l]  > $min+$binsize*$j){
	    $j++;
	}
	$hist[$i][$j]++;
	if (!defined $rep[$i][$j]){
	    $rep[$i][$j]=$l;
	}
	else{
	    my $d1=$xi[$i][$l]-($min-$binsize*$j);
	    my $d2=$rep[$i][$j]-($min+$binsize*$j);
	    if ($d1< $d2){ #$xi[$i][$l] is closer to the middle of the bin
		$rep[$i][$j]=$l;
	    }
	}
    }
#    for ($j=0; $j<$bins; $j++){
#	print "$i $j $hist[$i][$j]\n";
#    }
}

#Initialize the array fact
for ($k=0; $k<=$Nwind; $k++){
    $fact[$k]=$nt[$k];
    $ebf[$k]=1;
}

#Iterate
for (my $n=1; $n<=$Niter; $n++){
    print "$n\n";
    for ($k=0; $k<=$Nwind; $k++){
	my $ebfk=0.0;
	for ($i=0; $i<=$Nwind; $i++){
	    #for ($l=0; $l<$nt[$i]; $l++){ #This is correct
	    for ($b=0; $b<$bins; $b++){
		my $bottom=0.0;
		for ($j=0; $j<=$Nwind; $j++){
		    #$bottom+=$ebw[$i][$l][$j]*$fact[$j];
		    if (defined $rep[$i][$b]){
			$bottom+=$ebw[$i][$rep[$i][$b]][$j]*$fact[$j];
		    }
		}
		#$ebfk+=$ebw[$i][$l][$k]/$bottom;
		if($hist[$i][$b]){
		    #$ebfk+=$hist[$i][$b]*$ebw[$i][$l][$k]/$bottom;
		    $ebfk+=$hist[$i][$b]*$ebw[$i][$rep[$i][$b]][$k]/$bottom;
		}
	    }
	}
	$ebf2[$k]=$ebfk;
	$ebf[$k]=1/($ebf[0]*$ebfk);
	$fact[$k]=$nt[$k]*$ebf[$k];
    }
    my $conv=1;
    for ($k=0; $k<=$Nwind; $k++){
	my $delta=abs($kBT*log($ebf[$k]/$ebf2[$k]));
	$conv=0 if ($delta>=$tol); 
	$ebf[$k]=$ebf2[0]/$ebf2[$k];
    }
    last if ($conv>0);
}


