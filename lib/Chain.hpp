//Sean M. Law

#ifndef CHAIN_H
#define CHAIN_H

#include "Residue.hpp"

class Chain {
  private:
    std::string id;
    Residue *start, *end;

  public:
    Chain();
};

#endif
