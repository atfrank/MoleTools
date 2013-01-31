//Sean M. Law

#include "Chain.hpp"

#include <vector>
using namespace std;

class Molecule {
  private:
    vector<Chain> chn;
    vector<Residue> res;
    vector<Atom> atm;

  public:
  
    int readPDB (string *ifile, int model=0);
};
