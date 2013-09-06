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
  std::cerr << "         [-pdb pdbFile mode1[:mode2[:...[:modeN]]]]" << std::endl; 
  exit(0);
}

int main (int argc, char **argv){

  int i;
  unsigned int j, k, l;
  std::string covar;
  std::string currArg;
  Analyze *anin;
  std::vector<unsigned int> mvec;
  std::vector<unsigned int> mval;
  std::vector<unsigned int> mode; 
  std::string pdb;
  unsigned int nrow, ncol;
  Molecule* avgmol;
  Molecule* eigmol;
  double eigval;
  Eigen::MatrixXd eigvec;
  Atom* a;
  unsigned int nmodel;

  covar.clear();
  mvec.clear();
  mval.clear();
  mode.clear();
  pdb.clear();
  avgmol=NULL;
  eigmol=NULL;
  nmodel=1;


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
    anin->preAnalysis(); //Reads covariance matrix and diagonalizes it
    //std::cout << anin->getEigen().eigenvalues() << std::endl;
    //std::cout << anin->getEigen().eigenvectors().col(anin->getEigen().eigenvalues().rows()-1) << std::endl;

    //Extract PDBs
    if (mode.size() > 0 && pdb.length() > 0){
      nrow=anin->getEigen().eigenvalues().rows();
      avgmol=Molecule::readPDB(pdb);
      std::cout << std::endl << "MODEL " << nmodel << std::endl;
      avgmol->writePDB();
      std::cout << "ENDMDL" << std::endl;
      nmodel++;
      for (j=0; j< mode.size(); j++){
        if (mode.at(j) > nrow){
          std::cerr << "Warning: Skipping unknown Mode " << mval.at(j) << std::endl;
        }
        else{
          eigmol=avgmol->clone(true,true);
          eigval=anin->getEigen().eigenvalues()[nrow-mode.at(j)];
          eigvec=anin->getEigen().eigenvectors().col(nrow-mode.at(j))*eigval;
          k=0;
          for (l=0; l< eigmol->getAtmVecSize(); l++){
            a=eigmol->getAtom(l);
            a->setX(a->getX()+eigvec(k));
            k++;
            a->setY(a->getY()+eigvec(k));
            k++;
            a->setZ(a->getZ()+eigvec(k));
            k++;
          }
          std::cout << std::endl << "MODEL " << nmodel << std::endl;
          eigmol->writePDB();
          std::cout << "ENDMDL" << std::endl;
          nmodel++;
          if (k != 3*eigmol->getAtmVecSize()){
            std::cerr << "Warning: The number of atoms in the PDB file (" << eigmol->getAtmVecSize();
            std::cerr << ") and covariance matrix " << k << " do not match!" << std::endl;
          }
          if (eigmol != NULL){
            delete eigmol;
            eigmol=NULL;
          }
        }
      }
      if (avgmol != NULL){
        delete avgmol;
      }
    }

    //Extract eigenvalues
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

    //Extract eigenvectors
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
