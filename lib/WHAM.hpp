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
    std::string cmd;
    std::vector<double> Fguess;
    std::vector<double> F;
    std::vector< std::vector< std::vector<double> > > expBVE; //Is dynamic and can be jagged
		std::vector< std::vector< std::vector<double> > > expBVxEx; //Is dynamic and can be jagged
    unsigned int nWindow;
    std::string fMeta;
    std::ifstream metaFile;
    std::istream *metainp;
    std::vector<unsigned int> bins;
    double tol;
    unsigned int maxIter;
    std::vector<double> T; //Beta(i) = E.size() 
    double targetT; //Target temp
    bool factorFlag;
    double factor; //For use with Molecular Transfer Model (MTM)
    std::vector< std::vector<std::string> > inps;

  public:
    WHAM ();
    void appendCmd(const std::string &str);
    void genWHAMInput();
    void readMetadata();
    void processVtot(); //Total biasing potentials
    void processEtot(); //Total potential (T-Rex)
    void processCoor(); //Reaction coord
    bool iterateWHAM();
    void fixTemp();

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
    void setNWindow(const unsigned int &nwin);
    void setNWindow(const int &nwin);

    std::string getMeta();
    unsigned int getTempSize();
    double getTemp(const int &element);
    std::string getCmd();
    unsigned int getNWindow();
};

#endif
