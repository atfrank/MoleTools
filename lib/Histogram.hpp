//Sean M. Law

#ifndef HISTOGRAM_H
#define HISTOGRAM_H

#include <vector>
#include <iostream>

class Histogram {
  private:
    std::vector<unsigned int> histo;
    std::vector<double> MAX; //Global Max
    std::vector<double> MIN; //Global Min
    std::vector< std::vector< double > > data;

  public:
    Histogram(const unsigned &ninpin, const unsigned int &ndimin);
    void updateMAXMIN(const std::vector<double> &sin);
    void appendData(const std::vector<double> &sin, const unsigned int &dimin);
    void genHistogram();
};

#endif
