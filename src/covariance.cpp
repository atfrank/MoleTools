//Sean M. Law

#include "Misc.hpp"

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <vector>

void usage(){
  std::cerr << "Usage:   Covariance [-options] <input(s)>" << std::endl;
  std::cerr << "Options: [-col col1[:col2[:..[:colN]]] || -col col1=colN=incr]]]" << std::endl;
  exit(0);
}

int main (int argc, char **argv){

  int i;
  unsigned int j, k;
  std::string currArg;
  std::vector<std::string> ifiles;
  std::ifstream inpFile;
  std::istream* finp;
  std::string line;
  std::vector<double> s;
  std::vector<unsigned int> cols;
  std::vector<unsigned int> range;
  unsigned int incr;
  std::vector<double> avg;
  std::vector<double> dx;
  unsigned int n;

  ifiles.clear();
  avg.clear();
  dx.clear();
  n=0;

  for (i=1; i<argc; i++){
    currArg=argv[i];
    if (currArg.compare("-h") == 0 || currArg.compare("-help") == 0){
      usage();
    }
    else if (currArg.compare("-col") == 0 || currArg.compare("-cols") == 0){
      currArg=argv[++i];
      if (currArg.find(":") != std::string::npos || Misc::isdigit(currArg)){
        Misc::splitNum(currArg, ":", cols, false);
      }
      else if (currArg.find("=") != std::string::npos){
        Misc::splitNum(currArg, "=", range, false);
        incr=1;
        if (range.size() >= 2){
          if (range.size() >= 3){
            incr=range.at(2);
          }
        }
        for (j=range.at(0); j<= range.at(1); j=j+incr){
          cols.push_back(j-1);
        }
      }
      else{
        std::cerr << std::endl << "Error: Unrecognized column specification format";
        std::cerr << std::endl << std::endl;
        usage();
      } 
    }
    else{
      ifiles.push_back(currArg);
    }
  }

  avg.resize(cols.size());
  dx.resize(cols.size());
  for (j=0; j< cols.size(); j++){
    avg.at(j)=0.0;
    dx.at(j)=0.0;
  }

  //Get average
  for (j=0; j< ifiles.size(); j++){
    inpFile.open(ifiles.at(j).c_str(), std::ios::in);
    finp=&inpFile;

    while(finp->good() && !(finp->eof())){
      getline(*finp, line);
      Misc::splitNum(line, " \t", s, false);
      for (k=0; k< cols.size(); k++){
        if (cols.at(k) < s.size()){
          avg.at(k)+=s.at(cols.at(k));
        }
        else{
          //Shouldn't be here
        }
      }
      n++; 
    }

    if (inpFile.is_open()){
      inpFile.close();
    }
  }

  //Get difference from average
  for (j=0; j< ifiles.size(); j++){
    inpFile.open(ifiles.at(j).c_str(), std::ios::in);
    finp=&inpFile;

    while(finp->good() && !(finp->eof())){
      getline(*finp, line);
      Misc::splitNum(line, " \t", s, false);
      for (k=0; k< cols.size(); k++){
        if (cols.at(k) < s.size()){
          dx.at(k)+=(s.at(cols.at(k))-avg.at(k));
        }
        else{
          //Shouldn't be here
        } 
      }
    }

    if (inpFile.is_open()){
      inpFile.close();
    }
  } 

  //Get average difference and construct covariance

  return 0;
}
