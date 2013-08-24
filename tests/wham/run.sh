#!/bin/csh

foreach i (`seq 0 1 15`)
  awk '{if (NR%10==0){print $0}}' data.t$i > t.$i
end
