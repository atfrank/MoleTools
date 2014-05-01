#Sean M. Law

include ../makefile.defs

LIBOBJS = Analyze.o Atom.o Bin.o BTree.o Chain.o DTree.o Histogram.o LARMORD.o LinAlg.o Misc.o Mol2.o Molecule.o PCASSO.o PDB.o Prmtop.o Residue.o Select.o Trajectory.o Coor.o WHAM.o

all:	libmoletools.$(LIBEXT)

libmoletools.$(LIBEXT):	$(LIBOBJS)
	ar rc libmoletools.$(LIBEXT) $(LIBOBJS)

clean:
	rm -rf *.o
	rm -rf *.a
	rm -rf *.so
	rm -rf *.dll
	rm -rf *stackdump
	rm -rf *~

%.o:	%.cpp
	$(CXX) $(CXXFLAGS) $(OPTFLAGS) $(DEBUG) $(INCFLAGS) -c $< -o $@
