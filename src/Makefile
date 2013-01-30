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
  $(LIBDIR)/Vector.hpp $(LIBDIR)/Atom.hpp $(LIBDIR)/PDB.hpp \
  $(OBJDIR)/Vector.o $(OBJDIR)/Atom.o $(OBJDIR)/PDB.o
	$(CC) $(CFLAGS) $(OPTFLAGS) -o $@ manipPDB.cpp \
    $(OBJDIR)/Vector.o $(OBJDIR)/Atom.o $(OBJDIR)/PDB.o

$(LIBDIR)/Vector.cpp: 

$(OBJDIR)/Vector.o: $(LIBDIR)/Vector.cpp $(LIBDIR)/Vector.hpp
	$(CC) $(CFLAGS) $(OPTFLAGS) -o $@ -c $(LIBDIR)/Vector.cpp

$(OBJDIR)/Atom.o: $(LIBDIR)/Atom.cpp $(LIBDIR)/Atom.hpp
	$(CC) $(CFLAGS) $(OPTFLAGS) -o $@ -c $(LIBDIR)/Atom.cpp

$(OBJDIR)/PDB.o: $(LIBDIR)/PDB.cpp $(LIBDIR)/PDB.hpp
	$(CC) $(CFLAGS) $(OPTFLAGS) -o $@ -c $(LIBDIR)/PDB.cpp

clean:
	rm -rf $(BINDIR)/*
	rm -rf $(OBJDIR)/*.o
	rm -rf $(LIBDIR)/*~
	rm -rf $(TSTDIR)/*~
	rm -rf *~
