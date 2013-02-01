//Sean M. Law
#ifndef PDB_H
#define PDB_H

#include "Molecule.hpp"

//#include <string>

class PDB {
  private:
    //All instances are declared as arguments to the functions below

  public:
    static int writePDB (Molecule& mol, std::string *ofile);
};

#endif
