//Sean M. Law

#ifndef LINALG_H
#define LINALG_H

#include <vector>
#include <algorithm>
#include <iostream>
#include <cmath>

class LinAlg {
  private:
    std::vector< std::vector<double> > U;
    std::vector< std::vector<double> > V;
    std::vector<double> s;
    
  public:
    void SVD (const std::vector< std::vector <double> > &M);
    void transpose ();
};

#endif
