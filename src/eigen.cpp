//Sean M. Law

#include "Molecule.hpp"
#include "Analyze.hpp"
#include "Misc.hpp"
#include "Trajectory.hpp"

#include <iostream>
#include <fstream>
#include <cstdlib>

void usage(){
  std::cerr << "Usage:   eigen [-options] <covarFILE>" << std::endl;
  std::cerr << "Options: [-vector mode1[:mode2[...[:modeN]]] | -vector modestart=modestop[=modeincr]]" << std::endl;
  std::cerr << "         [-value mode1[:mode2[...[:modeN]]] | -value modestart=modestop[=modeincr]]" << std::endl;
  std::cerr << "         [-pdb PDBfile mode1[:mode2[:...[:modeN]]]] [-sel selection]" << std::endl;
  std::cerr << "         [-out TRAJfile]" << std::endl; 
  std::cerr << "         [-overlap covarFile]" << std::endl;
  exit(0);
}

int main (int argc, char **argv){

  int i;
  unsigned int j;
  std::string covar;
  std::string cmpcovar;
  std::string currArg;
  Analyze *anin;
  Analyze *ancmp;
  std::vector<unsigned int> mvec;
  std::vector<unsigned int> mval;
  std::vector<unsigned int> mode; 
  std::vector<unsigned int> range;
  unsigned int incr;
  std::string pdb;
  unsigned int nrow, ncol;
  Molecule* avgmol;
  unsigned int nmodel;
  std::string sel;
  std::string fout;
  std::ofstream trjout;
  Trajectory *ftrjout;

  covar.clear();
  cmpcovar.clear();
  anin=NULL;
  ancmp=NULL;
  mvec.clear();
  mval.clear();
  mode.clear();
  pdb.clear();
  avgmol=NULL;
  nmodel=0;
  sel=":.";
  ftrjout=NULL;


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
    }
    else if (currArg.compare("-sel") == 0){
      currArg=argv[++i];
      sel=currArg;
    }
    else if (currArg.compare("-out") == 0){
      currArg=argv[++i];
      fout=currArg;
    }
    else if (currArg.compare("-overlap") == 0){
      currArg=argv[++i];
      cmpcovar=currArg;
    }
    else if (currArg.compare("-vector") == 0){
      currArg=argv[++i];
      if (currArg.find(":") != std::string::npos || Misc::isdigit(currArg)){
        Misc::splitNum(currArg, ":", mvec, false);
      }
      else if (currArg.find("=") != std::string::npos){
        Misc::splitNum(currArg, "=", range, false);
        incr=1;
        if (range.size() >= 2){
          if (range.size() >= 3){
            incr=range.at(2);
          }
          for (j=range.at(0); j<= range.at(1); j=j+incr){
            mvec.push_back(j);
          }
        }
      }
      else{
        std::cerr << std::endl << "Error: Unrecongized mode format";
        std::cerr << std::endl << std::endl;
        usage();
      }
    }
    else if (currArg.compare("-value") == 0){
      currArg=argv[++i];
      if (currArg.find(":") != std::string::npos || Misc::isdigit(currArg)){
        Misc::splitNum(currArg, ":", mval, false);
      }
      else if (currArg.find("=") != std::string::npos){
        Misc::splitNum(currArg, "=", range, false);
        incr=1;
        if (range.size() >= 2){
          if (range.size() >= 3){
            incr=range.at(2);
          }
          for (j=range.at(0); j<= range.at(1); j=j+incr){
            mval.push_back(j);
          }
        }
      }
      else{
        std::cerr << std::endl << "Error: Unrecongized mode format";
        std::cerr << std::endl << std::endl;
        usage();
      }
    }
    else{
      covar=currArg;
    }
  }

  if (covar.length() > 0){
    anin=new AnalyzeCovariance;
    anin->setInput(covar);

    if (mode.size() > 0 && pdb.length() > 0){
      //Extract mode structures
      anin->addSel(sel);
      anin->preAnalysis(Molecule::readPDB(pdb)); //Reads covariance matrix and diagonalizes it
      nrow=anin->getEigen().eigenvalues().rows();
      if (fout.length() > 0){
        trjout.open(fout.c_str(), std::ios::binary);
        ftrjout=new Trajectory;
        ftrjout->setDefaultHeader();
        ftrjout->setMolecule(anin->getMol(0));
        ftrjout->setNAtom(static_cast<int>(anin->getMol(0)->getNAtomSelected()));
        ftrjout->writeHeader(trjout);
        ftrjout->writeFrame(trjout); //write input PDB
      } 
      else{
        std::cout << std::endl << "MODEL " << nmodel+1 << std::endl;
        anin->getMol(0)->writePDB(); //write input PDB
        std::cout << "ENDMDL" << std::endl;
      }
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
          std::cerr << "Warning: Skipping unknown Mode " << mode.at(j) << std::endl;
        }
        else{
          anin->setEigenMode(mode.at(j));
          if (trjout.is_open()){
            ftrjout->setMolecule(anin->getMol(1));
            ftrjout->writeFrame(trjout);
          }
          else{
            std::cout << std::endl << "MODEL " << nmodel+1 << std::endl;
            anin->getMol(1)->writePDB();
            std::cout << "ENDMDL" << std::endl;
          }
          nmodel++;
        }
      }
      if (avgmol != NULL){
        delete avgmol;
      }
      if (trjout.is_open()){
        ftrjout->setNFrame(nmodel);
        ftrjout->setNStep(nmodel);
        ftrjout->writeHeader(trjout);
        trjout.close();
        if (ftrjout != NULL){
          delete ftrjout;
        }
      }
    }
    else if (cmpcovar.length() > 0){
      //Eigenvector Overlap
      anin->preAnalysis();
      ancmp=new AnalyzeCovariance;
      ancmp->setInput(cmpcovar);
      ancmp->preAnalysis();
      anin->writeEigenOverlap(ancmp, mvec);
    }
    else if (mval.size() > 0){
      //Extract Eigenvalues
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
    else if (mvec.size() > 0){
      //Extract eigenvectors
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
    else{
      //Do nothing
    }
  }

  return 0;
}
