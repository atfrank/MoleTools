//Sean M. Law

#include "Analyze.hpp"
#include "Misc.hpp"

#include <iostream>
#include <fstream>
#include <cstdlib>

void usage(){
  std::cerr << "Usage:   eigen [-options] <covarFILE>" << std::endl;
  std::cerr << "Options: [-vector mode1[:mode2[...[:modeN]]] | -value mode1[:mode2[...[:modeN]]]]" << std::endl;
  exit(0);
}

int main (int argc, char **argv){

  int i;
  unsigned int j;
  std::string covar;
  std::string currArg;
  Analyze *anin;
  std::vector<unsigned int> mvec;
  std::vector<unsigned int> mval;
  unsigned int nrow, ncol;

  covar.clear();
  mvec.clear();
  mval.clear();


  for (i=1; i<argc; i++){
    currArg=argv[i];
    if (currArg.compare("-h") == 0 || currArg.compare("-help") == 0){
      usage();
    }
    else if (currArg.compare("-vector") == 0){
      currArg=argv[++i];
      Misc::splitNum(currArg, ":", mvec, false);
      mval.clear();
    }
    else if (currArg.compare("-value") == 0){
      currArg=argv[++i];
      Misc::splitNum(currArg, ":", mval, false);
      mvec.clear();
    }
    else{
      covar=currArg;
    }
  }

  if (covar.length() > 0){
    anin=new AnalyzeCovariance;
    anin->setInput(covar);
    anin->preAnalysis(); //Reads covariance matrix and diagonalizes it
    //std::cout << anin->getEigen().eigenvalues() << std::endl;
    //std::cout << anin->getEigen().eigenvectors().col(anin->getEigen().eigenvalues().rows()-1) << std::endl;
    if (mval.size() > 0){
      if (mval.size() == 1 && mval.at(0) == 0){
        //If mode = 0 then print all eigenvalues
        mval.clear();
        for (j=0; j< static_cast<unsigned int>(anin->getEigen().eigenvalues().rows()); j++){
          mval.push_back(j+1);
        }
      }
      nrow=anin->getEigen().eigenvalues().rows();
      for (j=0; j< mval.size(); j++){
        if (mval.at(j) > nrow){
          std::cerr << "Warning: Skipping unknown Mode " << mval.at(j) << std::endl;
        }
        else{
          std::cout << anin->getEigen().eigenvalues().row(nrow-mval.at(j)) << std::endl;
        }
      }
    }
    if (mvec.size() > 0){
      if (mvec.size() == 1 && mvec.at(0) == 0){
        //If mode = 0 then print all eigenvectors
        mvec.clear();
        for (j=0; j< static_cast<unsigned int>(anin->getEigen().eigenvectors().cols()); j++){
          mvec.push_back(j+1);
        }
      }
      ncol=anin->getEigen().eigenvectors().cols();
      for (j=0; j< mvec.size(); j++){
        if (mvec.at(j) > ncol){
          std::cerr << "Warning: Skipping unknown Mode " << mvec.at(j) << std::endl;
        }
        else{
          std::cout << anin->getEigen().eigenvectors().col(ncol-mvec.at(j)) << std::endl << std::endl;
        }
      }
    }
  }

  return 0;
}
