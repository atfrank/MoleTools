//Sean M. Law
#ifndef PDB_H
#define PDB_H

#include "Molecule.hpp"
#include "Misc.hpp"

#include <iostream>
#include <iomanip>
#include <map>
#include <stdexcept>

class PDB {
  private:
    std::map<std::string, int> chnMap;
		std::string format; //Output format

  public:
    PDB();
    static void writePDBFormat (Molecule* mol, std::ostringstream &out, bool selFlag=true, bool chnFlag=false);
    static Molecule* readPDB (std::string ifile, int model=0);
    Atom* processAtomLine (std::string line, Atom* lastAtom);
};

#endif
