#Sean M. Law

LIBDIR = ../lib
BINDIR = ../bin
TSTDIR = ../tests
OBJDIR = obj

CXX = g++
CXXFLAGS = -Wall -lm -I$(LIBDIR) 
#DEBUG = -g
OPTFLAGS = -O3

all: $(BINDIR) $(OBJDIR) modPDB modTraj secondary contact pdb2Traj traj2PDB;

$(BINDIR):
	mkdir -p $(BINDIR)

$(OBJDIR):
	mkdir -p $(OBJDIR)

secondary: $(BINDIR)/secondary;
$(BINDIR)/secondary: $(OBJDIR)/secondary.o $(OBJDIR)/SABA.o \
		$(OBJDIR)/Molecule.o \
		$(OBJDIR)/Chain.o $(OBJDIR)/Residue.o \
		$(OBJDIR)/Atom.o $(OBJDIR)/Vector.o \
		$(OBJDIR)/Constants.o $(OBJDIR)/PDB.o \
		$(OBJDIR)/Select.o $(OBJDIR)/Misc.o \
		$(OBJDIR)/Analyze.o
	$(CXX) $(CXXFLAGS) $(OPTFLAGS) $(DEBUG) -o $@ \
		$(OBJDIR)/secondary.o $(OBJDIR)/SABA.o \
		$(OBJDIR)/Molecule.o $(OBJDIR)/Chain.o \
		$(OBJDIR)/Residue.o $(OBJDIR)/Atom.o \
		$(OBJDIR)/Vector.o $(OBJDIR)/Constants.o \
		$(OBJDIR)/PDB.o $(OBJDIR)/Select.o \
		$(OBJDIR)/Misc.o $(OBJDIR)/Analyze.o

contact: $(BINDIR)/contact;
$(BINDIR)/contact: $(OBJDIR)/contact.o $(OBJDIR)/Trajectory.o \
		$(OBJDIR)/Molecule.o $(OBJDIR)/Chain.o \
		$(OBJDIR)/Residue.o $(OBJDIR)/Atom.o \
		$(OBJDIR)/Vector.o $(OBJDIR)/Constants.o \
		$(OBJDIR)/PDB.o $(OBJDIR)/Select.o \
		$(OBJDIR)/Misc.o $(OBJDIR)/Analyze.o
	$(CXX) $(CXXFLAGS) $(OPTFLAGS) $(DEBUG) -o $@ \
		$(OBJDIR)/contact.o $(OBJDIR)/Trajectory.o \
		$(OBJDIR)/Molecule.o $(OBJDIR)/Chain.o \
		$(OBJDIR)/Residue.o $(OBJDIR)/Atom.o \
		$(OBJDIR)/Vector.o $(OBJDIR)/Constants.o \
		$(OBJDIR)/PDB.o $(OBJDIR)/Select.o \
		$(OBJDIR)/Misc.o $(OBJDIR)/Analyze.o

pdb2Traj: $(BINDIR)/pdb2Traj;
$(BINDIR)/pdb2Traj: $(OBJDIR)/pdb2Traj.o $(OBJDIR)/Trajectory.o \
		$(OBJDIR)/Molecule.o $(OBJDIR)/Chain.o \
		$(OBJDIR)/Residue.o $(OBJDIR)/Atom.o \
		$(OBJDIR)/Vector.o $(OBJDIR)/Constants.o \
		$(OBJDIR)/PDB.o $(OBJDIR)/Select.o \
		$(OBJDIR)/Misc.o $(OBJDIR)/Analyze.o
	$(CXX) $(CXXFLAGS) $(OPTFLAGS) $(DEBUG) -o $@ \
		$(OBJDIR)/pdb2Traj.o $(OBJDIR)/Trajectory.o \
		$(OBJDIR)/Molecule.o $(OBJDIR)/Chain.o \
		$(OBJDIR)/Residue.o $(OBJDIR)/Atom.o \
		$(OBJDIR)/Vector.o $(OBJDIR)/Constants.o \
		$(OBJDIR)/PDB.o $(OBJDIR)/Select.o \
		$(OBJDIR)/Misc.o $(OBJDIR)/Analyze.o

traj2PDB: $(BINDIR)/traj2PDB;
$(BINDIR)/traj2PDB: $(OBJDIR)/traj2PDB.o $(OBJDIR)/Trajectory.o \
		$(OBJDIR)/Molecule.o $(OBJDIR)/Chain.o \
		$(OBJDIR)/Residue.o $(OBJDIR)/Atom.o \
		$(OBJDIR)/Vector.o $(OBJDIR)/Constants.o \
		$(OBJDIR)/PDB.o $(OBJDIR)/Select.o \
		$(OBJDIR)/Misc.o $(OBJDIR)/Analyze.o
	$(CXX) $(CXXFLAGS) $(OPTFLAGS) $(DEBUG) -o $@ \
		$(OBJDIR)/traj2PDB.o $(OBJDIR)/Trajectory.o \
		$(OBJDIR)/Molecule.o $(OBJDIR)/Chain.o \
		$(OBJDIR)/Residue.o $(OBJDIR)/Atom.o \
		$(OBJDIR)/Vector.o $(OBJDIR)/Constants.o \
		$(OBJDIR)/PDB.o $(OBJDIR)/Select.o \
		$(OBJDIR)/Misc.o $(OBJDIR)/Analyze.o

modTraj: $(BINDIR)/modTraj;
$(BINDIR)/modTraj: $(OBJDIR)/modTraj.o $(OBJDIR)/Trajectory.o \
		$(OBJDIR)/Molecule.o $(OBJDIR)/Chain.o \
		$(OBJDIR)/Residue.o $(OBJDIR)/Atom.o \
		$(OBJDIR)/Vector.o $(OBJDIR)/Constants.o \
		$(OBJDIR)/PDB.o $(OBJDIR)/Select.o \
		$(OBJDIR)/Misc.o $(OBJDIR)/Analyze.o
	$(CXX) $(CXXFLAGS) $(OPTFLAGS) $(DEBUG) -o $@ \
		$(OBJDIR)/modTraj.o $(OBJDIR)/Trajectory.o \
		$(OBJDIR)/Molecule.o $(OBJDIR)/Chain.o \
		$(OBJDIR)/Residue.o $(OBJDIR)/Atom.o \
		$(OBJDIR)/Vector.o $(OBJDIR)/Constants.o \
		$(OBJDIR)/PDB.o $(OBJDIR)/Select.o \
		$(OBJDIR)/Misc.o $(OBJDIR)/Analyze.o

modPDB: $(BINDIR)/modPDB;
$(BINDIR)/modPDB: $(OBJDIR)/modPDB.o $(OBJDIR)/Molecule.o \
		$(OBJDIR)/Chain.o $(OBJDIR)/Residue.o \
		$(OBJDIR)/Atom.o $(OBJDIR)/Vector.o \
		$(OBJDIR)/Constants.o $(OBJDIR)/PDB.o \
		$(OBJDIR)/Select.o $(OBJDIR)/Misc.o \
		$(OBJDIR)/Analyze.o
	$(CXX) $(CXXFLAGS) $(OPTFLAGS) $(DEBUG) -o $@ \
		$(OBJDIR)/modPDB.o $(OBJDIR)/Molecule.o \
		$(OBJDIR)/Chain.o $(OBJDIR)/Residue.o \
		$(OBJDIR)/Atom.o $(OBJDIR)/Vector.o \
		$(OBJDIR)/Constants.o $(OBJDIR)/PDB.o \
		$(OBJDIR)/Select.o $(OBJDIR)/Misc.o \
		$(OBJDIR)/Analyze.o

$(OBJDIR)/secondary.o: secondary.cpp \
		$(LIBDIR)/Molecule.hpp $(LIBDIR)/SABA.hpp
	$(CXX) $(CXXFLAGS) $(OPTFLAGS) $(DEBUG) -o $@ -c secondary.cpp

$(OBJDIR)/traj2PDB.o: traj2PDB.cpp \
		$(LIBDIR)/Molecule.hpp $(LIBDIR)/Trajectory.hpp
	$(CXX) $(CXXFLAGS) $(OPTFLAGS) $(DEBUG) -o $@ -c traj2PDB.cpp

$(OBJDIR)/pdb2Traj.o: pdb2Traj.cpp \
		$(LIBDIR)/Molecule.hpp $(LIBDIR)/Trajectory.hpp
	$(CXX) $(CXXFLAGS) $(OPTFLAGS) $(DEBUG) -o $@ -c pdb2Traj.cpp

$(OBJDIR)/contact.o: contact.cpp \
		$(LIBDIR)/Molecule.hpp $(LIBDIR)/Trajectory.hpp $(LIBDIR)/Misc.hpp \
		$(LIBDIR)/Analyze.hpp
	$(CXX) $(CXXFLAGS) $(OPTFLAGS) $(DEBUG) -o $@ -c contact.cpp

$(OBJDIR)/modTraj.o: modTraj.cpp \
		$(LIBDIR)/Molecule.hpp $(LIBDIR)/Trajectory.hpp
	$(CXX) $(CXXFLAGS) $(OPTFLAGS) $(DEBUG) -o $@ -c modTraj.cpp

$(OBJDIR)/modPDB.o: modPDB.cpp \
		$(LIBDIR)/Molecule.hpp $(LIBDIR)/Vector.hpp 
	$(CXX) $(CXXFLAGS) $(OPTFLAGS) $(DEBUG) -o $@ -c modPDB.cpp

$(OBJDIR)/Trajectory.o: $(LIBDIR)/Trajectory.cpp $(LIBDIR)/Trajectory.hpp \
		$(LIBDIR)/Molecule.hpp
	$(CXX) $(CXXFLAGS) $(OPTFLAGS) $(DEBUG) -o $@ -c $(LIBDIR)/Trajectory.cpp

$(OBJDIR)/PDB.o: $(LIBDIR)/PDB.cpp $(LIBDIR)/PDB.hpp \
		$(LIBDIR)/Molecule.hpp 
	$(CXX) $(CXXFLAGS) $(OPTFLAGS) $(DEBUG) -o $@ -c $(LIBDIR)/PDB.cpp

$(OBJDIR)/Analyze.o: $(LIBDIR)/Analyze.cpp $(LIBDIR)/Analyze.hpp \
		$(LIBDIR)/Molecule.hpp $(LIBDIR)/Vector.hpp $(LIBDIR)/Constants.hpp
	$(CXX) $(CXXFLAGS) $(OPTFLAGS) $(DEBUG) -o $@ -c $(LIBDIR)/Analyze.cpp

$(OBJDIR)/Molecule.o: $(LIBDIR)/Molecule.cpp $(LIBDIR)/Molecule.hpp \
		$(LIBDIR)/PDB.hpp $(LIBDIR)/Chain.hpp $(LIBDIR)/Analyze.hpp 
	$(CXX) $(CXXFLAGS) $(OPTFLAGS) $(DEBUG) -o $@ -c $(LIBDIR)/Molecule.cpp

$(OBJDIR)/Chain.o: $(LIBDIR)/Chain.cpp $(LIBDIR)/Chain.hpp \
		$(LIBDIR)/Residue.hpp 
	$(CXX) $(CXXFLAGS) $(OPTFLAGS) $(DEBUG) -o $@ -c $(LIBDIR)/Chain.cpp

$(OBJDIR)/Residue.o: $(LIBDIR)/Residue.cpp $(LIBDIR)/Residue.hpp \
		$(LIBDIR)/Atom.hpp 
	$(CXX) $(CXXFLAGS) $(OPTFLAGS) $(DEBUG) -o $@ -c $(LIBDIR)/Residue.cpp

$(OBJDIR)/Atom.o: $(LIBDIR)/Atom.cpp $(LIBDIR)/Atom.hpp \
		$(LIBDIR)/Vector.hpp 
	$(CXX) $(CXXFLAGS) $(OPTFLAGS) $(DEBUG) -o $@ -c $(LIBDIR)/Atom.cpp

$(OBJDIR)/Vector.o: $(LIBDIR)/Vector.cpp $(LIBDIR)/Vector.hpp 
	$(CXX) $(CXXFLAGS) $(OPTFLAGS) $(DEBUG) -o $@ -c $(LIBDIR)/Vector.cpp

$(OBJDIR)/Constants.o: $(LIBDIR)/Constants.cpp $(LIBDIR)/Constants.hpp
	$(CXX) $(CXXFLAGS) $(OPTFLAGS) $(DEBUG) -o $@ -c $(LIBDIR)/Constants.cpp

$(OBJDIR)/Select.o: $(LIBDIR)/Select.cpp $(LIBDIR)/Select.hpp \
		$(LIBDIR)/Molecule.hpp
	$(CXX) $(CXXFLAGS) $(OPTFLAGS) $(DEBUG) -o $@ -c $(LIBDIR)/Select.cpp

$(OBJDIR)/Misc.o: $(LIBDIR)/Misc.cpp $(LIBDIR)/Misc.hpp
	$(CXX) $(CXXFLAGS) $(OPTFLAGS) $(DEBUG) -o $@ -c $(LIBDIR)/Misc.cpp

$(OBJDIR)/SABA.o: $(LIBDIR)/SABA.cpp $(LIBDIR)/SABA.hpp
	$(CXX) $(CXXFLAGS) $(OPTFLAGS) $(DEBUG) -o $@ -c $(LIBDIR)/SABA.cpp

clean:
	rm -rf $(BINDIR)/*
	rm -rf $(OBJDIR)/*.o
	rm -rf $(LIBDIR)/*~
	rm -rf $(TSTDIR)/*~
	rm -rf *~
	rm -rf *stackdump
