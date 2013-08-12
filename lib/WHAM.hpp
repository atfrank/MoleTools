// Sean M. Law

#ifndef WHAM_H
#define WHAM_H

#include "Misc.hpp"
#include "Constants.hpp"

#include <vector>

class WHAM {
  private:
    std::vector<double> Fguess;
    std::vector<double> F;
    std::vector< std::vector< std::vector<double> > > E; //Is dynamic and can be jagged
		std::vector< std::vector< std::vector<double> > > Ex; //Is dynamic and can be jagged
    std::string fMeta;
    std::vector<unsigned int> bins;
    double tol;
    unsigned int maxIter;

  public:
    WHAM ();
    void genWHAMInput();
    void readWHAMInput();
    void iterateWHAM();

    void setMeta(const std::string &metain);
    void setBins(const std::string &binsin);
    void setBins(const std::vector<unsigned int> &binsin);
    void setBins(const std::vector<int> &binsin);
    void setTol(const double &tolin);
    void setMaxIter(const unsigned int &iterin);
};

#endif
