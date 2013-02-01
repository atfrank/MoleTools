#Sean M. Law

LIBDIR = ../lib
BINDIR = ../bin
TSTDIR = ../test
OBJDIR = obj

CC = g++
CFLAGS = -Wall -lm -I$(LIBDIR)
OPTFLAGS = -O2

all: $(BINDIR) $(OBJDIR) manipPDB;

$(BINDIR):
	mkdir -p $(BINDIR)

$(OBJDIR):
	mkdir -p $(OBJDIR)

manipPDB: $(BINDIR)/manipPDB;
$(BINDIR)/manipPDB: manipPDB.cpp \
		$(LIBDIR)/PDB.hpp $(OBJDIR)/PDB.o
		$(CC) $(CFLAGS) $(OPTFLAGS) -o $@ manipPDB.cpp \
			$(OBJDIR)/Molecule.o $(OBJDIR)/Chain.o $(OBJDIR)/Residue.o \
      $(OBJDIR)/Atom.o $(OBJDIR)/Vector.o $(OBJDIR)/PDB.o

$(OBJDIR)/PDB.o: $(LIBDIR)/PDB.cpp $(LIBDIR)/PDB.hpp \
		$(LIBDIR)/Molecule.hpp $(OBJDIR)/Molecule.o
		$(CC) $(CFLAGS) $(OPTFLAGS) -o $@ -c $(LIBDIR)/PDB.cpp

$(OBJDIR)/Molecule.o: $(LIBDIR)/Molecule.cpp $(LIBDIR)/Molecule.hpp \
		$(LIBDIR)/PDB.hpp $(LIBDIR)/Chain.hpp \
		$(OBJDIR)/PDB.o $(OBJDIR)/Chain.o
		$(CC) $(CFLAGS) $(OPTFLAGS) -o $@ -c $(LIBDIR)/Molecule.cpp

$(OBJDIR)/Chain.o: $(LIBDIR)/Chain.cpp $(LIBDIR)/Chain.hpp \
		$(LIBDIR)/Residue.hpp $(OBJDIR)/Residue.o
		$(CC) $(CFLAGS) $(OPTFLAGS) -o $@ -c $(LIBDIR)/Chain.cpp

$(OBJDIR)/Residue.o: $(LIBDIR)/Residue.cpp $(LIBDIR)/Residue.hpp \
		$(LIBDIR)/Atom.hpp $(OBJDIR)/Atom.o
		$(CC) $(CFLAGS) $(OPTFLAGS) -o $@ -c $(LIBDIR)/Residue.cpp

$(OBJDIR)/Atom.o: $(LIBDIR)/Atom.cpp $(LIBDIR)/Atom.hpp \
		$(LIBDIR)/Vector.hpp $(OBJDIR)/Vector.o
		$(CC) $(CFLAGS) $(OPTFLAGS) -o $@ -c $(LIBDIR)/Atom.cpp

$(OBJDIR)/Vector.o: $(LIBDIR)/Vector.cpp $(LIBDIR)/Vector.hpp
	$(CC) $(CFLAGS) $(OPTFLAGS) -o $@ -c $(LIBDIR)/Vector.cpp

clean:
	rm -rf $(BINDIR)/*
	rm -rf $(OBJDIR)/*.o
	rm -rf $(LIBDIR)/*~
	rm -rf $(TSTDIR)/*~
	rm -rf *~
