//Sean M. Law

#ifndef RESIDUE_H
#define RESIDUE_H

#include "Atom.hpp"

#include <vector>

class Residue {
  private:
    std::string *resname;
    int *resid;
    std::string *chainid;
    Atom *start;
    Atom *end;
    std::string *segid;
    std::vector<Atom*> atmVec;
    bool sel;

  public:
    Residue();

    void reset();
    int getResId();
    std::string getResName();
    std::string getChainId();
    Atom* getStart();
    Atom* getEnd();
    std::string getSegId();
    void addAtom(Atom* atmEntry);
    Atom* getAtom (int element);
    unsigned int getAtmVecSize();
    void setSel(bool selin);
    bool& getSel();
    void selAll();
    void deselAll();
};

#endif
