#!/bin/sh

#../bin/manipPDB -nsel A:ARG._B:140-145.CA ../tests/2DHB.pdb  > t.pdb
#../bin/manipPDB -nsel !A+B:1-137._!A+B:135+ARG. ../tests/2DHB.pdb  > t.pdb
#../bin/manipPDB -nsel A:.CA ../tests/2DHB.pdb  > t.pdb
#../bin/manipPDB -nsel :^100-135. ../tests/2DHB.pdb  > t.pdb
../bin/manipPDB -nsel A+B+\ :. ../tests/2DHB.pdb > t.pdb

