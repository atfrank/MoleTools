//Sean M. Law

#include "WHAM.hpp"

#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>

#define MAXINPFILES 4096
#define MAXLINESIZE 4096

void usage(){
  std::cerr << "Usage:   wham [-options] <metadatafile>" << std::endl;
  std::cerr << "Options: [-bins rcoor1[:rcoor2[:rcoor3]]]" << std::endl;
  std::cerr << "         [-iter value] [-tol value | -Ftol value]" << std::endl;
  std::cerr << "         [-temp value]" << std::endl;
  std::cerr << std::endl;
  exit(0);
}

int main (int argc, char **argv){

  int i;
  std::string currArg;
  std::vector<unsigned int> bins;
  double tol;
  unsigned int maxIter;
  WHAM *wham;
 
  tol=1E-5;
  maxIter=1E6;
  wham=new WHAM;

  for (i=1; i<argc; i++){
    currArg=argv[i];
    if (currArg == "-h" || currArg == "-help"){
      usage();
    }
    else if (currArg == "-bins"){
      currArg=argv[++i];
      wham->setBins(currArg);
    }
    else if (currArg == "-iter" || currArg == "-maxiter"){
      currArg=argv[++i];
      std::stringstream(currArg) >> maxIter;
      wham->setMaxIter(maxIter);
    }
    else if (currArg == "-tol" || currArg == "-ftol"){
      currArg=argv[++i];
      std::stringstream(currArg) >> tol;
      wham->setTol(tol);
    }
    else if (currArg == "-Ftol"){
      currArg=argv[++i];
      std::stringstream(currArg) >> tol;
      tol=log(tol);
      wham->setTol(tol);
    }
    else{
      wham->setMeta(currArg);
    }
  }


  return 0;
}
