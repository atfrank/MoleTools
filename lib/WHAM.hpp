// Sean M. Law

#ifndef WHAM_H
#define WHAM_H

#include "Misc.hpp"
#include "Constants.hpp"

#include <iostream>
#include <fstream>
#include <vector>

class WHAM {
  private:
    std::vector<double> Fguess;
    std::vector<double> F;
    std::vector< std::vector< std::vector<double> > > E; //Is dynamic and can be jagged
		std::vector< std::vector< std::vector<double> > > Ex; //Is dynamic and can be jagged
    std::string fMeta;
    std::ifstream metaFile;
    std::istream *metainp;
    std::vector<unsigned int> bins;
    double tol;
    unsigned int maxIter;
    std::vector<double> T; //Beta(i) + Beta(0) = E.size() + 1
    bool factorFlag;
    double factor; //For use with Molecular Transfer Model (MTM)

  public:
    WHAM ();
    void genWHAMInput();
    void processMeta();
    void iterateWHAM();

    void setMeta(const std::string &metain);
    void setBins(const std::string &binsin);
    void setBins(const std::vector<unsigned int> &binsin);
    void setBins(const std::vector<int> &binsin);
    void setTol(const double &tolin);
    void setMaxIter(const unsigned int &iterin);
    bool setTemp(const std::string &tin);
    void setTemp(const std::vector<double> &tin);
    bool setTempRange(const std::string &tin);
    void setFactor(const double &factorin);

    std::string getMeta();
    unsigned int getTempSize();
    double getTemp(const int &element);
};

#endif
