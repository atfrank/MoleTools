//Sean M. Law

#ifndef CHAIN_H
#define CHAIN_H

#include "Residue.hpp"

#include <vector>

class Chain {
  private:
    std::string id;
    std::vector<Atom> atmVec;

  public:
    Chain();

    void reset();

    //Get Chain info
    int getChainIdSize();
    std::string getChainId();

    //Set Chain info
    void setChainId(const std::string& chainidin);
    void setChainId(); //Clear
};

#endif
