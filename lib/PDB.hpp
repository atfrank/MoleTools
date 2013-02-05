//Sean M. Law
#ifndef PDB_H
#define PDB_H

#include "Molecule.hpp"

#include <map>

class PDB {
  private:
    std::string lastChain;
    bool newChain;
    int lastRes;
    bool ter;
    std::map<std::string, int> chnMap;

  public:
    PDB();
    void setLastChain(std::string chainid);
    void setNewChain(bool val);
    void setLastRes(int resid);
    static std::string writePDBFormat (Molecule& mol);
    static void readPDB (Molecule& mol, std::string ifile, int model=0);
    Atom processAtomLine (std::string line);
    void processResidue (Molecule& mol);
    void processChain (Molecule& mol);
    void setTer (bool val);
};

#endif
