// Sean M. Law

#ifndef WHAM_H
#define WHAM_H

#include "Misc.hpp"
#include "Constants.hpp"

#include <vector>

class WHAM {
  private:
    unsigned int nSims;
    std::vector<unsigned int> nData;
    std::vector<double> Fguess;
    std::vector< std::vector< std::vector<double> > > energy; //Is dynamic and can be jagged
    unsigned int bins;
    double tol;
    unsigned int maxIter;
    bool extrapFlag; //
    bool guessFlag; //

  public:
    WHAM ();
    void genWHAMInput();
    void readWHAMInput();
    void iterateWHAM();
};

#endif
