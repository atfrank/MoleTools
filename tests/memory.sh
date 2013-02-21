#!/bin/csh

#make clean
#make

#valgrind --tool=memcheck --dsymutil=yes --leak-check=full --log-file=memory.log ../bin/manipPDB ../tests/1SB0.pdb

#valgrind --tool=memcheck --dsymutil=yes --leak-check=full --show-reachable=yes --log-file=memory.log ../bin/manipPDB ../tests/1SB0.pdb

valgrind --tool=memcheck --leak-check=full --show-reachable=yes --log-file=memory.log ../bin/manipPDB ../tests/1SB0.pdb
