//Sean M. Law

#ifndef CHAIN_H
#define CHAIN_H

#include "Residue.hpp"

#include <vector>

class Chain {
  private:
    std::string *id;
    std::vector<Residue *> resVec; //Vector of residue pointers
    std::vector<Atom *> atmVec; //Vector of atom pointers
//    bool sel;

  public:
    Chain();

    void reset();

    void addResidue(Residue* resEntry);
    void addAtom(Atom* atmEntry);

    //Get Chain info
    Atom* getAtom(int element);
    Residue* getResidue(int element);
    std::string getChainId();
    unsigned int getAtmVecSize();
    unsigned int getResVecSize();
//    void setSel(bool selin);
//    bool& getSel();
    void selAll();
    void deselAll();
};

#endif
