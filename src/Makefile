#Sean M. Law

include ../makefile.defs

OBJS = libtools.$(LIBEXT)

all:	$(BINDIR) libtools.$(LIBEXT) unitTest contact modPDB modTraj pdb2Traj secondary traj2PDB correl analyzeTraj

clean:
	rm -rf *.o
	rm -rf *stackdump
	rm -rf ../bin/*

$(BINDIR):
	mkdir -p $(BINDIR)

libtools.$(LIBEXT):
	$(MAKE) -C ../lib libtools.$(LIBEXT)

unitTest: unitTest.o $(OBJS)
	$(CXX) $(CXXFLAGS) $(OPTFLAGS) $(DEBUG) $(LINKFLAGS) -o $(BINDIR)/$@ $<

contact: contact.o $(OBJS)
	$(CXX) $(CXXFLAGS) $(OPTFLAGS) $(DEBUG) $(LINKFLAGS) -o $(BINDIR)/$@ $<

modPDB: modPDB.o $(OBJS)
	$(CXX) $(CXXFLAGS) $(OPTFLAGS) $(DEBUG) $(LINKFLAGS) -o $(BINDIR)/$@ $<

modTraj: modTraj.o $(OBJS)
	$(CXX) $(CXXFLAGS) $(OPTFLAGS) $(DEBUG) $(LINKFLAGS) -o $(BINDIR)/$@ $<

pdb2Traj: pdb2Traj.o $(OBJS)
	$(CXX) $(CXXFLAGS) $(OPTFLAGS) $(DEBUG) $(LINKFLAGS) -o $(BINDIR)/$@ $<

secondary: secondary.o $(OBJS)
	$(CXX) $(CXXFLAGS) $(OPTFLAGS) $(DEBUG) $(LINKFLAGS) -o $(BINDIR)/$@ $<

traj2PDB: traj2PDB.o $(OBJS)
	$(CXX) $(CXXFLAGS) $(OPTFLAGS) $(DEBUG) $(LINKFLAGS) -o $(BINDIR)/$@ $<

correl: correl.o $(OBJS)
	$(CXX) $(CXXFLAGS) $(OPTFLAGS) $(DEBUG) $(LINKFLAGS) -o $(BINDIR)/$@ $<

analyzeTraj: analyzeTraj.o $(OBJS)
	$(CXX) $(CXXFLAGS) $(OPTFLAGS) $(DEBUG) $(LINKFLAGS) -o $(BINDIR)/$@ $<

%.o:	%.cpp
	$(CXX) $(CXXFLAGS) $(OPTFLAGS) $(DEBUG) -c $(INCFLAGS) $<  
