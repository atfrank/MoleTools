//Sean M. Law

#include "Chain.hpp"

#include <vector>
using namespace std;

class Molecule {
  private:
    vector<Chain> chnVec;
    vector<Residue> resVec;
    vector<Atom> atmVec;

  public:
  
    int readPDB (string *ifile, int model=0);
    Atom getAtom(int atmnum);
};
