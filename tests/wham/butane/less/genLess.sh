#!/bin/csh

set k = 0.02

rm metadata

foreach i (`seq 1 1 72`)
  awk -v k=$k '{if (NR%10==0){for(j=1; j<=72; j++){x0=180-(j-1)*5; if (x0<0){x0=x0+360}; printf "%s ", k*($2-x0)*($2-x0)} printf "\n"}}' ../grossfield.$i > v.$i
  awk '{if (NR%10==0){print $2}}' ../grossfield.$i > rc.$i
  echo "v.$i 0 rc.$i" >> metadata
end
