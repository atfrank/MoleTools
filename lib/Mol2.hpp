//Sean M. Law
#ifndef MOL2_H
#define MOL2_H

#include "Molecule.hpp"
#include "Misc.hpp"
#include "Select.hpp"

#include <iostream>
#include <iomanip>
#include <map>
#include <stdexcept>

class Mol2 {
  private:
    std::map<std::string, int> chnMap;

  public:
    Mol2();
    static void writeMol2Format (Molecule* mol, std::ostringstream &out, bool selFlag=true);
    static Molecule* readMol2 (std::string ifile);
    Atom* processAtomLine (std::string line, Atom* lastAtom);
};

#endif
