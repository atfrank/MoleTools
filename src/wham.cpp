//Sean M. Law

#include "WHAM.hpp"

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cstdlib>

#define MAXINPFILES 4096
#define MAXLINESIZE 4096

void usage(){
  std::cerr << "Usage:   wham [-options] <metadatafile>" << std::endl;
  std::cerr << "Options: [-bins rcoor1[:rcoor2[:rcoor3]]]" << std::endl;
  std::cerr << "         [-iter value] [-tol value | -Ftol value]" << std::endl;
  std::cerr << "         [-temp T1[:T2...[[:TN[:Ttarget]]]] | -temp T1=TN=[incr[=Ttarget]]]" << std::endl;
  std::cerr << "         [-fguess file | -fval file]" << std::endl;
  std::cerr << "         [-factor value]" << std::endl;
  std::cerr << std::endl;
  exit(0);
}

int main (int argc, char **argv){

  int i;
  unsigned int j;
  std::string currArg;
  std::vector<unsigned int> bins;
  double tol;
  unsigned int maxIter;
  double factor;
  WHAM *wham;
  std::string fguess;
  std::string fval;
 
  tol=1E-5;
  maxIter=1E6;
  factor=1.0;
  wham=new WHAM;
  j=0;
  fguess.clear();
  fval.clear();

  //Copy original command
  for (i=0; i<argc; i++){
    wham->appendCmd(std::string(" ")+argv[i]);
  }

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
    else if (currArg == "-temp" || currArg == "-temps"){
      currArg=argv[++i];
      if (currArg.find(":") != std::string::npos || Misc::isdouble(currArg)){
        if (wham->setTemp(currArg) == true){
          usage();
        }
      }
      else if (currArg.find("=") != std::string::npos){
        if (wham->setTempRange(currArg) == true){
          usage();
        }
      }
      else{
        std::cerr << std::endl << "Error: Unrecognized temperature format";
        std::cerr << std::endl << std::endl;
        usage();
      }
    }
    else if (currArg == "-fval"){
      currArg=argv[++i];
      fval=currArg;
      fguess.clear();
    }
    else if (currArg == "-fguess"){
      currArg=argv[++i];
      fguess=currArg;
      fval.clear();
    }
    else if (currArg == "-factor"){
      currArg=argv[++i];
      std::stringstream(currArg) >> factor;
      wham->setFactor(factor);
    }
    else{
      wham->setMeta(currArg);
    }
  }

  if (wham->getMeta().length() <= 0){
    std::cerr << std::endl << "Error:  Please provide an input metadata file";
    std::cerr << std::endl << std::endl;
    usage();
  }

  wham->readMetadata();
  wham->processEnergies();
  if (fguess.size() != 0){
    wham->setFguess(fguess);
  }
  if (fval.size() != 0){
    wham->setFval(fval); //No WHAM iteration needed
  }
  else{
    wham->iterateWHAM();
    std::cout << "#" << wham->getCmd() << std::endl;
  }

  wham->setDenomInv();

  if (wham->processCoor() == true){ //Reaction coordinates
     wham->binOnTheFly();
     wham->printPMF();
  }
  else{
    std::cerr << "Error: The number of datapoints did not match up" << std::endl;
  }
  
//  for (j=0; j< wham->getTempSize(); j++){
//    std::cerr << std::fixed << std::setprecision(6) << wham->getTemp(j) << std::endl;
//  }

  return 0;
}
