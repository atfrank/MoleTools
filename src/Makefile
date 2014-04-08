#Sean M. Law

include ../makefile.defs

OBJS = libmoletools.$(LIBEXT)

all:	$(BINDIR) libmoletools.$(LIBEXT) analyzeTraj contact contactOrder correl covariance eigen histo mine modPDB modTraj mol22PDB pcasso pdb2Traj traj2PDB unitTest wham

clean:
	rm -rf *.o
	rm -rf *stackdump
	rm -rf ../bin/*
	rm -rf *~

$(BINDIR):
	mkdir -p $(BINDIR)

libmoletools.$(LIBEXT):
	$(MAKE) -C ../MoleTools libmoletools.$(LIBEXT)

analyzeTraj: analyzeTraj.o $(OBJS)
	$(CXX) $(CXXFLAGS) $(OPTFLAGS) $(DEBUG) $(LINKFLAGS) -o $(BINDIR)/$@ $<

contact: contact.o $(OBJS)
	$(CXX) $(CXXFLAGS) $(OPTFLAGS) $(DEBUG) $(LINKFLAGS) -o $(BINDIR)/$@ $<

contactOrder: contactOrder.o $(OBJS)
	$(CXX) $(CXXFLAGS) $(OPTFLAGS) $(DEBUG) $(LINKFLAGS) -o $(BINDIR)/$@ $<

correl: correl.o $(OBJS)
	$(CXX) $(CXXFLAGS) $(OPTFLAGS) $(DEBUG) $(LINKFLAGS) -o $(BINDIR)/$@ $<
	
covariance: covariance.o $(OBJS)
	$(CXX) $(CXXFLAGS) $(OPTFLAGS) $(DEBUG) $(LINKFLAGS) -o $(BINDIR)/$@ $<

eigen: eigen.o $(OBJS)
	$(CXX) $(CXXFLAGS) $(OPTFLAGS) $(DEBUG) $(LINKFLAGS) -o $(BINDIR)/$@ $<

histo: histo.o $(OBJS)
	$(CXX) $(CXXFLAGS) $(OPTFLAGS) $(DEBUG) $(LINKFLAGS) -o $(BINDIR)/$@ $<

mine: mine.o $(OBJS)
	$(CXX) $(CXXFLAGS) $(OPTFLAGS) $(DEBUG) -L../MINE -lmine $(LINKFLAGS) -o $(BINDIR)/$@ $<

modPDB: modPDB.o $(OBJS)
	$(CXX) $(CXXFLAGS) $(OPTFLAGS) $(DEBUG) $(LINKFLAGS) -o $(BINDIR)/$@ $<

modTraj: modTraj.o $(OBJS)
	$(CXX) $(CXXFLAGS) $(OPTFLAGS) $(DEBUG) $(LINKFLAGS) -o $(BINDIR)/$@ $<

mol22PDB: mol22PDB.o $(OBJS)
	$(CXX) $(CXXFLAGS) $(OPTFLAGS) $(DEBUG) $(LINKFLAGS) -o $(BINDIR)/$@ $<
	
pcasso: pcasso.o $(OBJS)
	$(CXX) $(CXXFLAGS) $(OPTFLAGS) $(DEBUG) $(LINKFLAGS) -o $(BINDIR)/$@ $<

pdb2Traj: pdb2Traj.o $(OBJS)
	$(CXX) $(CXXFLAGS) $(OPTFLAGS) $(DEBUG) $(LINKFLAGS) -o $(BINDIR)/$@ $<

secondary: secondary.o $(OBJS)
	$(CXX) $(CXXFLAGS) $(OPTFLAGS) $(DEBUG) $(LINKFLAGS) -o $(BINDIR)/$@ $<

traj2PDB: traj2PDB.o $(OBJS)
	$(CXX) $(CXXFLAGS) $(OPTFLAGS) $(DEBUG) $(LINKFLAGS) -o $(BINDIR)/$@ $<

unitTest: unitTest.o $(OBJS)
	$(CXX) $(CXXFLAGS) $(OPTFLAGS) $(DEBUG) $(LINKFLAGS) -o $(BINDIR)/$@ $<

wham: wham.o $(OBJS)
	$(CXX) $(CXXFLAGS) $(OPTFLAGS) $(DEBUG) $(LINKFLAGS) -o $(BINDIR)/$@ $<

%.o:	%.cpp
	$(CXX) $(CXXFLAGS) $(OPTFLAGS) $(DEBUG) -c $(INCFLAGS) $<
