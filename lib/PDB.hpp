//Sean M. Law
#ifndef PDB_H
#define PDB_H

#include "Molecule.hpp"

//#include <string>

class PDB {
  private:
    //All instances are declared as arguments to the functions below

  public:
    static std::string writePDBFormat (Molecule& mol);
    static void readPDB (Molecule& mol, std::string ifile, int model=0);
    static int getCurrModel (std::string line);
    static Atom processAtomLine (std::string line);
};

#endif
