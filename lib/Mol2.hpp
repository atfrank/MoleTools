//Sean M. Law
#ifndef MOL2_H
#define MOL2_H

#include <map>
#include <string>

class Molecule;
class Atom;

class Mol2 {
  private:
    std::map<std::string, int> chnMap;

  public:
    Mol2();
    static Molecule* readMol2 (const std::string ifile, const std::string format="");
    Atom* processAtomLine (const std::string line, Atom* lastAtom);
};

#endif
