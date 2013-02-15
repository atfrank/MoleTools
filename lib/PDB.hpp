//Sean M. Law
#ifndef PDB_H
#define PDB_H

#include "Molecule.hpp"

#include <map>
#include <iostream>
#include <iomanip>
#include <stdexcept>

class PDB {
  private:
    std::map<std::string, int> chnMap;

  public:
    PDB();
    static std::string writePDBFormat (Molecule* mol);
    static Molecule* readPDB (std::string ifile, int model=0);
    Atom* processAtomLine (std::string line, Atom* lastAtom);
};

#endif
