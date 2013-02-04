//Sean M. Law

#ifndef RESIDUE_H
#define RESIDUE_H

#include "Atom.hpp"

#include <vector>

class Residue {
  private:
    std::string resname;
    int resnum;
    std::string chainid;
    Atom *start;
    Atom *end;
    std::string segid;
    std::vector<Atom *> atmVec;

  public:
    Residue();
    void setResName(std::string resnamein);
    void setResNum(int resnumin);
    void setChainId(std::string chainidin);
    void setStart(Atom* startin);
    void setEnd(Atom* endin);
    void setSegId(std::string segidin);
    int getResNum();
    void addAtom(Atom* atmEntry);
};

#endif
