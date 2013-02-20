//Sean M. Law

#ifndef SABA_H
#define SABA_H

#include <Molecule.hpp>

class SABA {
  private:
    std::vector<std::string> ss;

  public:
		static Molecule* getPseudoCenter(Molecule *mol);
};

#endif
