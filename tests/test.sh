#!/bin/sh

#../bin/manipPDB -nsel A:ARG._B:140-145.CA ../tests/1E3M.pdb  > t.pdb
#../bin/manipPDB -nsel !A+B:1-137._!A+B:135+ARG. ../tests/1E3M.pdb  > t.pdb
#../bin/manipPDB -nsel A:.CA ../tests/1E3M.pdb  > t.pdb
#../bin/manipPDB -nsel :^100-135. ../tests/1E3M.pdb  > t.pdb
../bin/manipPDB -sel A:.CA_B:.CA ../tests/1E3M.pdb > t.pdb

