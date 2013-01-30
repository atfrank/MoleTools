//Sean M. Law

#include "Chain.hpp"

class Molecule {
  private:
    Chain *chn;
    Residue *res;
    Atom *atm;

  public:
  
    int readPDB (string *ifile, int model=0);
};
