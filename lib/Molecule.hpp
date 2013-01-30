//Sean M. Law

#include "Chain.hpp"

#include <vector>

class Molecule {
  private:
    vector<Chain> chn;
    vector<Residue> res;
    vector<Atom> atm;

  public:
  
    int readPDB (string *ifile, int model=0);
};
