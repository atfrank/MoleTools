#!/bin/sh

#../bin/modPDB -nsel A:ARG._B:140-145.CA ../tests/1E3M.pdb  > t.pdb
#../bin/modPDB -nsel !A+B:1-137._!A+B:135+ARG. ../tests/1E3M.pdb  > t.pdb
#../bin/modPDB -nsel A:.CA ../tests/1E3M.pdb  > t.pdb
#../bin/modPDB -nsel :^100-135. ../tests/1E3M.pdb  > t.pdb
../bin/modPDB -sel A:.CA_B:.CA ../tests/1E3M.pdb > t.pdb
bin/modPDB -data tests/data.txt -occup -atom tests/file.ref.pdb
bin/modPDB -data tests/data.txt -bfac -atom tests/file.ref.pdb
bin/modPDB -data tests/data.txt -atom -residue tests/file.ref.pdb
bin/modPDB -data tests/data.txt -occup -residue tests/file.ref.pdb
