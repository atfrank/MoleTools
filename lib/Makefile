#Sean M. Law

include ../makefile.defs

LIBOBJS = LinAlg.o Analyze.o Atom.o Chain.o Constants.o Misc.o Molecule.o PDB.o Residue.o SABA.o Select.o Trajectory.o Vector.o Mol2.o

all:	libtools.$(LIBEXT)

libtools.$(LIBEXT):	$(LIBOBJS)
	$(CXX) $(CXXFLAGS) $(OPTFLAGS) $(LIBFLAGS) $^ -o libtools.$(LIBEXT)

clean:
	rm -rf *.o
	rm -rf *.$(LIBEXT)
	rm -rf *stackdump

%.o:	%.cpp
	$(CXX) $(CXXFLAGS) $(OPTFLAGS) $(LIBFLAGS) $(DEBUG) -c $< -o $@
