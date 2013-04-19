#Sean M. Law

include ../makefile.defs

LIBOBJS = Analyze.o Atom.o Chain.o Constants.o Misc.o Molecule.o PDB.o Residue.o SABA.o Select.o Trajectory.o Vector.o

all:	libtools.so

libtools.so:	$(LIBOBJS)
	$(CXX) $(CXXFLAGS) $(OPTFLAGS) $(LIBFLAGS) $^ -o libtools.so

clean:
	rm -rf *.o
	rm -rf *.so
	rm -rf *stackdump

%.o:	%.cpp
	$(CXX) $(CXXFLAGS) $(OPTFLAGS) $(LIBFLAGS) $(DEBUG) -c $< -o $@
