//Sean M. Law
//Aaron T. Frank
    
/*
This file is part of MoleTools.

MoleTools is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MoleTools is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MoleTools.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "Molecule.hpp"
#include "Analyze.hpp"
#include "Misc.hpp"
#include "Trajectory.hpp"

#include <fstream>

void usage(){
  std::cerr << "Usage:   eigen [-options] <covarFILE>" << std::endl;
  std::cerr << "Options: [-vector] [-value]" << std::endl;
  std::cerr << "         [-mode mode1[:mode2[...[:modeN]]] | -mode modestart=modestop[=modeincr]]" << std::endl;
  std::cerr << "         [-pdb PDBfile] [-sel selection]" << std::endl;
  std::cerr << "         [-out TRAJfile]" << std::endl;
  std::cerr << "         [-entropy | -quasi] [-temp value]" << std::endl;
  std::cerr << "         [-top file]" << std::endl;
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
  bool vector;
  bool value;
  std::vector<unsigned int> modes;
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
  std::string top;
  bool quasi;
  bool entropy;
  double temp;

  covar.clear();
  cmpcovar.clear();
  anin=NULL;
  ancmp=NULL;
  vector=false;
  value=false;
  modes.clear();
  pdb.clear();
  avgmol=NULL;
  nmodel=0;
  sel=":.";
  ftrjout=NULL;
  top.clear();
  entropy=false;
  quasi=false;
  temp=300;


  for (i=1; i<argc; i++){
    currArg=argv[i];
    if (currArg.compare("-h") == 0 || currArg.compare("-help") == 0){
      usage();
    }
    else if (currArg.compare("-pdb") == 0){
      currArg=argv[++i];
      pdb=currArg;
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
    else if (currArg.compare("-vector") == 0 || currArg.compare("-vectors") == 0){
      vector=true;
    }
    else if (currArg.compare("-value") == 0 || currArg.compare("-values") == 0){
      value=true;
    }
    else if (currArg.compare("-mode") == 0 || currArg.compare("-modes") == 0){
      currArg=argv[++i];
      if (currArg.find(":") != std::string::npos || Misc::isdigit(currArg)){
        Misc::splitNum(currArg, ":", modes, false);
      }
      else if (currArg.find("=") != std::string::npos){
        Misc::splitNum(currArg, "=", range, false);
        incr=1;
        if (range.size() >= 2){
          if (range.size() >= 3){
            incr=range.at(2);
          }
          for (j=range.at(0); j<= range.at(1); j=j+incr){
            modes.push_back(j);
          }
        }
      }
      else{
        std::cerr << std::endl << "Error: Unrecongized mode format";
        std::cerr << std::endl << std::endl;
        usage();
      }
    }
    else if (currArg.compare("-entropy") == 0){
      entropy=true;
      quasi=false;
    }
    else if (currArg.compare("-quasi") == 0){
      quasi=true;
      entropy=false;
    }
    else if (currArg.compare("-temp") == 0){
      currArg=argv[++i];
      std::stringstream(currArg) >> temp;
    }
    else if (currArg.compare("-top") == 0){
      currArg=argv[++i];
      top=currArg;
    }
		else if (currArg.compare(0,1,"-") == 0){
      std::cerr << "Warning: Skipping unknown option \"" << currArg << "\"" << std::endl;
    }
    else{
      covar=currArg;
    }
  }

  if (quasi == true && pdb.length() == 0){
    std::cerr << "Error: Please provide a PDB file via \"-pdb\" for quasiharmonic entropy calculation" << std::endl;
    usage();
  }

  if (covar.length() > 0){
    anin=new AnalyzeCovariance;
    anin->setInput(covar);

    if (pdb.length() > 0){
      anin->addSel(sel);
      anin->preAnalysis(Molecule::readPDB(pdb), top); //Reads covariance matrix and diagonalizes it
      nrow=anin->getEigen().eigenvalues().rows();
      if (quasi == true){
        //Extract entropy
        if (modes.size() == 0 || (modes.size() == 1 && modes.at(0) == 0)){
          //If modes = 0 then print all models
          modes.clear();
          for (j=0; j< static_cast<unsigned int>(anin->getEigen().eigenvalues().rows()); j++){
            modes.push_back(j+1);
          }
        }
        std::cout << anin->quasiharmonicEntropy(anin->getMol(0), anin->getEigen(), modes, temp) << std::endl;
      }
      else{
        //Extract mode structures
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
        if (modes.size() == 0 || (modes.size() == 1 && modes.at(0) == 0)){
          //If modes = 0 then print all models
          modes.clear();
          for (j=0; j< static_cast<unsigned int>(anin->getEigen().eigenvalues().rows()); j++){
            modes.push_back(j+1);
          }
        }
        for (j=0; j< modes.size(); j++){
          if (modes.at(j) > nrow){
            std::cerr << "Warning: Skipping unknown Mode " << modes.at(j) << std::endl;
          }
          else{
            anin->setEigenMode(modes.at(j));
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
    }
    else if (entropy == true){
      anin->preAnalysis();
      if ((modes.size() == 1 && modes.at(0) == 0) || modes.size() == 0){
        //If mode = 0 then print all eigenvalues
        modes.clear();
        for (j=0; j< static_cast<unsigned int>(anin->getEigen().eigenvalues().rows()); j++){
          modes.push_back(j+1);
        }
      }
      std::cout << anin->configurationalEntropy(anin->getEigen(), modes) << std::endl;
    }
    else if (cmpcovar.length() > 0){
      //Eigenvector Overlap
      anin->preAnalysis();
      ancmp=new AnalyzeCovariance;
      ancmp->setInput(cmpcovar);
      ancmp->preAnalysis();
      anin->writeEigenOverlap(ancmp, modes);
    }
    else if (value == true){
      //Extract Eigenvalues
      anin->preAnalysis();
      if (modes.size() == 1 && modes.at(0) == 0){
        //If mode = 0 then print all eigenvalues
        modes.clear();
        for (j=0; j< static_cast<unsigned int>(anin->getEigen().eigenvalues().rows()); j++){
          modes.push_back(j+1);
        }
      }
      nrow=anin->getEigen().eigenvalues().rows();
      for (j=0; j< modes.size(); j++){
        if (modes.at(j) > nrow){
          std::cerr << "Warning: Skipping unknown Mode " << modes.at(j) << std::endl;
        }
        else{
          std::cout << anin->getEigen().eigenvalues().row(nrow-modes.at(j)) << std::endl;
        }
      }
    }
    else if (vector == true){
      //Extract eigenvectors
      anin->preAnalysis();
      if (modes.size() == 1 && modes.at(0) == 0){
        //If mode = 0 then print all eigenvectors
        modes.clear();
        for (j=0; j< static_cast<unsigned int>(anin->getEigen().eigenvectors().cols()); j++){
          modes.push_back(j+1);
        }
      }
      ncol=anin->getEigen().eigenvectors().cols();
      for (j=0; j< modes.size(); j++){
        if (modes.at(j) > ncol){
          std::cerr << "Warning: Skipping unknown Mode " << modes.at(j) << std::endl;
        }
        else{
          std::cout << anin->getEigen().eigenvectors().col(ncol-modes.at(j)) << std::endl << std::endl;
        }
      }
    }
    else{
      //Do nothing
    }
  }

  return 0;
}
