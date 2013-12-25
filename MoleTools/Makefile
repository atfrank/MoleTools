#Sean M. Law

include ../makefile.defs

LIBOBJS = Analyze.o Atom.o Chain.o Constants.o Enum.o Histogram.o LinAlg.o Misc.o Molecule.o PDB.o Prmtop.o Residue.o SABA.o Select.o BTree.o Trajectory.o Vector.o Mol2.o WHAM.o

all:	libmoletools.$(LIBEXT)

libmoletools.$(LIBEXT):	$(LIBOBJS)
	$(CXX) $(CXXFLAGS) $(OPTFLAGS) $(LIBFLAGS) $^ -o libmoletools.$(LIBEXT)

clean:
	rm -rf *.o
	rm -rf *.$(LIBEXT)
	rm -rf *stackdump
	rm -rf *~

%.o:	%.cpp
	$(CXX) $(CXXFLAGS) $(OPTFLAGS) $(LIBFLAGS) $(DEBUG) -c $< -o $@
