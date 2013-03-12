//Sean M. Law
#ifndef PDB_H
#define PDB_H

#include "Molecule.hpp"
#include "Misc.hpp"

#include <iostream>
#include <iomanip>
#include <map>
#include <stdexcept>
#include <iostream>

class PDB {
  private:
    std::map<std::string, int> chnMap;

  public:
    PDB();
    static std::string writePDBFormat (Molecule* mol, bool selFlag=true);
    static Molecule* readPDB (std::string ifile, int model=0);
    Atom* processAtomLine (std::string line, Atom* lastAtom);
};

#endif
