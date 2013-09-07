//Sean M. Law

#include "Molecule.hpp"
#include "Analyze.hpp"
#include "Misc.hpp"

#include <iostream>
#include <fstream>
#include <cstdlib>

void usage(){
  std::cerr << "Usage:   eigen [-options] <covarFILE>" << std::endl;
  std::cerr << "Options: [-vector mode1[:mode2[...[:modeN]]] | -value mode1[:mode2[...[:modeN]]]]" << std::endl;
  std::cerr << "         [-pdb pdbFile mode1[:mode2[:...[:modeN]]]] [-sel selection]" << std::endl; 
  std::cerr << "         [-comp covarFile]" << std::endl;
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
  std::vector<unsigned int> mode; 
  std::string pdb;
  unsigned int nrow, ncol;
  Molecule* avgmol;
  unsigned int nmodel;
  std::string sel;

  covar.clear();
  mvec.clear();
  mval.clear();
  mode.clear();
  pdb.clear();
  avgmol=NULL;
  nmodel=1;
  sel=":.";


  for (i=1; i<argc; i++){
    currArg=argv[i];
    if (currArg.compare("-h") == 0 || currArg.compare("-help") == 0){
      usage();
    }
    else if (currArg.compare("-pdb") == 0){
      currArg=argv[++i];
      pdb=currArg;
      currArg=argv[++i];
      Misc::splitNum(currArg, ":", mode, false);
      mval.clear();
      mvec.clear();
    }
    else if (currArg.compare("-sel") == 0){
      currArg=argv[++i];
      sel=currArg;
    }
    else if (currArg.compare("-vector") == 0){
      currArg=argv[++i];
      Misc::splitNum(currArg, ":", mvec, false);
      mval.clear();
      mode.clear();
    }
    else if (currArg.compare("-value") == 0){
      currArg=argv[++i];
      Misc::splitNum(currArg, ":", mval, false);
      mvec.clear();
      mode.clear();
    }
    else{
      covar=currArg;
    }
  }

  if (covar.length() > 0){
    anin=new AnalyzeCovariance;
    anin->setInput(covar);

    //Extract PDBs
    if (mode.size() > 0 && pdb.length() > 0){
      anin->addSel(sel);
      anin->preAnalysis(Molecule::readPDB(pdb)); //Reads covariance matrix and diagonalizes it
      nrow=anin->getEigen().eigenvalues().rows();
      std::cout << std::endl << "MODEL " << nmodel << std::endl;
      anin->getMol(0)->writePDB(); //write input PDB
      std::cout << "ENDMDL" << std::endl;
      nmodel++;
      if (mode.size() == 1 && mode.at(0) == 0){
        //If mode = 0 then print all models
        mode.clear();
        for (j=0; j< static_cast<unsigned int>(anin->getEigen().eigenvalues().rows()); j++){
          mode.push_back(j+1);
        }
      }
      for (j=0; j< mode.size(); j++){
        if (mode.at(j) > nrow){
          std::cerr << "Warning: Skipping unknown Mode " << mval.at(j) << std::endl;
        }
        else{
          std::cout << std::endl << "MODEL " << nmodel << std::endl;
          anin->setEigenMode(mode.at(j));
          anin->getMol(1)->writePDB();
          std::cout << "ENDMDL" << std::endl;
          nmodel++;
        }
      }
      if (avgmol != NULL){
        delete avgmol;
      }
    }

    //Extract eigenvalues
    if (mval.size() > 0){
      anin->preAnalysis();
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

    //Extract eigenvectors
    if (mvec.size() > 0){
      anin->preAnalysis();
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
