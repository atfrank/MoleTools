//Sean M. Law

#include "Histogram.hpp"

Histogram::Histogram(const unsigned &ninpin, const unsigned int &ndimin){
  MIN.resize(ndimin);
  MAX.resize(ndimin);
  for (unsigned int i=0; i< ndimin; i++){
    MIN.at(i)=1E20;
    MAX.at(i)=-1E20;
  }
  data.resize(ninpin);
}

void Histogram::updateMAXMIN(const std::vector<double> &sin){
  //Global Max/Min
  for (unsigned int i=0; i< MIN.size(); i++){
    if (sin.at(i) < MIN.at(i)){
      MIN.at(i)=sin.at(i);
    }
    if (sin.at(i) > MAX.at(i)){
      MAX.at(i)=sin.at(i);
    }
  }
}

void Histogram::appendData(const std::vector<double> &sin, const unsigned int &dimin){
  this->updateMAXMIN(sin);
  data.at(dimin).push_back(sin);
//  data.push_back(sin);
}

void Histogram::genHistogram(){
//  for (unsigned int i=0; i< data
}

