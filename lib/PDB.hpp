//Sean M. Law
#ifndef PDB_H
#define PDB_H

#include "Molecule.hpp"

#include <map>

class PDB {
  private:
    std::string lastChain;
    int ter;
    std::map<std::string, int> chnMap;

  public:
    PDB();
    static std::string writePDBFormat (Molecule& mol);
    static void readPDB (Molecule& mol, std::string ifile, int model=0);
    Atom processAtomLine (std::string line);
    Chain processChain(Atom atmEntry);
    void setTer (int terin);
};

#endif
