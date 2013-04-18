include ../makefile.defs

OBJS = libtools.so

all:	$(BINDIR) libtools.so contact modPDB modTraj pdb2Traj secondary traj2PDB

clean:
	rm -rf *.o
	rm -rf *stackdump

$(BINDIR):
	mkdir -p $(BINDIR)

libtools.so:
	$(MAKE) -C ../lib libtools.so

contact: contact.o $(OBJS)
	$(CXX) $(CXXFLAGS) $(DEBUG) $(LINKFLAGS) -o $(BINDIR)/$@ $<

modPDB: modPDB.o $(OBJS)
	$(CXX) $(CXXFLAGS) $(DEBUG) $(LINKFLAGS) -o $(BINDIR)/$@ $<

modTraj: modTraj.o $(OBJS)
	$(CXX) $(CXXFLAGS) $(DEBUG) $(LINKFLAGS) -o $(BINDIR)/$@ $<

pdb2Traj: pdb2Traj.o $(OBJS)
	$(CXX) $(CXXFLAGS) $(DEBUG) $(LINKFLAGS) -o $(BINDIR)/$@ $<

secondary: secondary.o $(OBJS)
	$(CXX) $(CXXFLAGS) $(DEBUG) $(LINKFLAGS) -o $(BINDIR)/$@ $<

traj2PDB: traj2PDB.o $(OBJS)
	$(CXX) $(CXXFLAGS) $(DEBUG) $(LINKFLAGS) -o $(BINDIR)/$@ $<

%.o:	%.cpp
	$(CXX) $(CXXFLAGS) $(DEBUG) -c $(INCFLAGS) $<  
