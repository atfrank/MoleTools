//Sean M. Law

#include "Chain.hpp"

class Molecule {
  private:
    Chain *chn;
    Residue *res;
    Atom *atm;

  int readPDB (string *ifile);
};
