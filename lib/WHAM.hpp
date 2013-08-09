// Sean M. Law

#ifndef WHAM_H
#define WHAM_H

#include "Misc.hpp"
#include "Constants.hpp"

#include <vector>

class WHAM {
  private:
    std::vector<double> Fguess;
    std::vector< std::vector< std::vector<double> > > E; //Is dynamic and can be jagged
		std::vector< std::vector< std::vector<double> > > Ex; //Is dynamic and can be jagged
    unsigned int bins;
    double tol;
    unsigned int maxIter;

  public:
    WHAM ();
    void genWHAMInput();
    void readWHAMInput();
    void iterateWHAM();
};

#endif
