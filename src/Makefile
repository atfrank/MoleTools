#Sean M. Law

LIBDIR = ../lib
BINDIR = ../bin
TSTDIR = ../test
OBJDIR = obj

CC = g++
CFLAGS = -Wall -lm -I$(LIBDIR)
#DEBUG = -g
OPTFLAGS = -O2

all: $(BINDIR) $(OBJDIR) manipPDB secondary;

$(BINDIR):
	mkdir -p $(BINDIR)

$(OBJDIR):
	mkdir -p $(OBJDIR)

secondary: $(BINDIR)/secondary;
$(BINDIR)/secondary: $(OBJDIR)/secondary.o $(OBJDIR)/SABA.o \
		$(OBJDIR)/Molecule.o \
		$(OBJDIR)/Chain.o $(OBJDIR)/Residue.o \
		$(OBJDIR)/Atom.o $(OBJDIR)/Vector.o \
		$(OBJDIR)/PDB.o $(OBJDIR)/Select.o \
		$(OBJDIR)/Misc.o
	$(CC) $(CFLAGS) $(OPTFLAGS) $(DEBUG) -o $@ \
		$(OBJDIR)/secondary.o $(OBJDIR)/SABA.o $(OBJDIR)/Molecule.o \
		$(OBJDIR)/Chain.o $(OBJDIR)/Residue.o $(OBJDIR)/Atom.o \
		$(OBJDIR)/Vector.o $(OBJDIR)/PDB.o $(OBJDIR)/Select.o \
		$(OBJDIR)/Misc.o


manipPDB: $(BINDIR)/manipPDB;
$(BINDIR)/manipPDB: $(OBJDIR)/manipPDB.o $(OBJDIR)/Molecule.o \
		$(OBJDIR)/Chain.o $(OBJDIR)/Residue.o \
		$(OBJDIR)/Atom.o $(OBJDIR)/Vector.o \
		$(OBJDIR)/PDB.o $(OBJDIR)/Select.o \
		$(OBJDIR)/Misc.o
	$(CC) $(CFLAGS) $(OPTFLAGS) $(DEBUG) -o $@ \
		$(OBJDIR)/manipPDB.o $(OBJDIR)/Molecule.o $(OBJDIR)/Chain.o \
		$(OBJDIR)/Residue.o $(OBJDIR)/Atom.o $(OBJDIR)/Vector.o \
		$(OBJDIR)/PDB.o $(OBJDIR)/Select.o $(OBJDIR)/Misc.o

$(OBJDIR)/secondary.o: secondary.cpp \
		$(LIBDIR)/Molecule.hpp $(LIBDIR)/SABA.hpp
	$(CC) $(CFLAGS) $(OPTFLAGS) $(DEBUG) -o $@ -c secondary.cpp

$(OBJDIR)/manipPDB.o: manipPDB.cpp \
		$(LIBDIR)/Molecule.hpp
	$(CC) $(CFLAGS) $(OPTFLAGS) $(DEBUG) -o $@ -c manipPDB.cpp

$(OBJDIR)/PDB.o: $(LIBDIR)/PDB.cpp $(LIBDIR)/PDB.hpp \
		$(LIBDIR)/Molecule.hpp 
	$(CC) $(CFLAGS) $(OPTFLAGS) $(DEBUG) -o $@ -c $(LIBDIR)/PDB.cpp

$(OBJDIR)/Molecule.o: $(LIBDIR)/Molecule.cpp $(LIBDIR)/Molecule.hpp \
		$(LIBDIR)/PDB.hpp $(LIBDIR)/Chain.hpp 
	$(CC) $(CFLAGS) $(OPTFLAGS) $(DEBUG) -o $@ -c $(LIBDIR)/Molecule.cpp

$(OBJDIR)/Chain.o: $(LIBDIR)/Chain.cpp $(LIBDIR)/Chain.hpp \
		$(LIBDIR)/Residue.hpp 
	$(CC) $(CFLAGS) $(OPTFLAGS) $(DEBUG) -o $@ -c $(LIBDIR)/Chain.cpp

$(OBJDIR)/Residue.o: $(LIBDIR)/Residue.cpp $(LIBDIR)/Residue.hpp \
		$(LIBDIR)/Atom.hpp 
	$(CC) $(CFLAGS) $(OPTFLAGS) $(DEBUG) -o $@ -c $(LIBDIR)/Residue.cpp

$(OBJDIR)/Atom.o: $(LIBDIR)/Atom.cpp $(LIBDIR)/Atom.hpp \
		$(LIBDIR)/Vector.hpp 
	$(CC) $(CFLAGS) $(OPTFLAGS) $(DEBUG) -o $@ -c $(LIBDIR)/Atom.cpp

$(OBJDIR)/Vector.o: $(LIBDIR)/Vector.cpp $(LIBDIR)/Vector.hpp
	$(CC) $(CFLAGS) $(OPTFLAGS) $(DEBUG) -o $@ -c $(LIBDIR)/Vector.cpp

$(OBJDIR)/Select.o: $(LIBDIR)/Select.cpp $(LIBDIR)/Select.hpp \
		$(LIBDIR)/Molecule.hpp
	$(CC) $(CFLAGS) $(OPTFLAGS) $(DEBUG) -o $@ -c $(LIBDIR)/Select.cpp

$(OBJDIR)/Misc.o: $(LIBDIR)/Misc.cpp $(LIBDIR)/Misc.hpp
	$(CC) $(CFLAGS) $(OPTFLAGS) $(DEBUG) -o $@ -c $(LIBDIR)/Misc.cpp

$(OBJDIR)/SABA.o: $(LIBDIR)/SABA.cpp $(LIBDIR)/SABA.hpp
	$(CC) $(CFLAGS) $(OPTFLAGS) $(DEBUG) -o $@ -c $(LIBDIR)/SABA.cpp

clean:
	rm -rf $(BINDIR)/*
	rm -rf $(OBJDIR)/*.o
	rm -rf $(LIBDIR)/*~
	rm -rf $(TSTDIR)/*~
	rm -rf *~
	rm -rf *stackdump
