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

#include "Trajectory.hpp"
#include "Molecule.hpp"

#include <iostream>
#include <fstream>

void usage(){
  std::cerr << std::endl << std::endl;
  std::cerr << "Usage:   pdb2Traj [-options] <PDBfile(s)>" << std::endl;
  std::cerr << "         [-out TRAJname] [-outsel selection]" << std::endl;
  exit(0);
}

int main (int argc, char **argv){

  int i;
  unsigned int j;
  std::string fout;
  Molecule *mol;
  std::vector<std::string> pdbs;
  std::string currArg;
  std::string outsel="";
  std::ofstream trjout;
  std::ifstream trjin;
  Trajectory *ftrjout;
  Trajectory *ftrjin;

  pdbs.clear();
  mol=NULL;
  ftrjout=NULL;
  ftrjin=NULL;

  for (i=1; i<argc; i++){
    currArg=argv[i];
    if (currArg.compare("-h") == 0 || currArg.compare("-help") == 0){
      usage();
    }
    else if (currArg.compare("-out") == 0){
      currArg=argv[++i];
      fout=currArg;
    }
    else if (currArg.compare("-outsel") == 0){
      currArg=argv[++i];
      outsel=currArg;
    }
    else if (currArg.compare(0,1,"-") == 0){
      std::cerr << "Warning: Skipping unknown option \"" << currArg << "\"" << std::endl;
    }
    else{
      pdbs.push_back(currArg);
    }
  }

  if (pdbs.size() == 0){
    std::cerr << std::endl << "Error: Please provide an input trajectory file" << std::endl << std::endl;
    usage();
  }
  else if (fout.length() == 0){
    std::cerr << std::endl << "Error: Please specify an output trajectory via \"-out\"" << std::endl; 
  }
  else{
    trjout.open(fout.c_str(), std::ios::binary);
    ftrjout=new Trajectory;
    ftrjout->setDefaultHeader();
    ftrjout->setNFrame(static_cast<int>(pdbs.size()));
    ftrjout->setNStep(static_cast<int>(pdbs.size()));

    for (j=0; j< pdbs.size(); j++){
      mol=Molecule::readPDB(pdbs.at(j));
      if (outsel.length() > 0){
        mol->select(outsel);
        mol=mol->clone(true, false);
      }
      
      if (j==0){
        ftrjout->setNAtom(static_cast<int>(mol->getNAtom()));
        ftrjout->writeHeader(trjout);
      }

      if (ftrjout->getNAtom() == static_cast<int>(mol->getNAtom())){
        ftrjout->setMolecule(mol);
  
        if (trjout.is_open()){
          if (trjin.is_open()){
            ftrjout->writeFrame(trjout, ftrjin);
          }
          else{
            ftrjout->writeFrame(trjout);
          }
        }
      }
      else{
        ftrjout->setNFrame(ftrjout->getNFrame()-1);
        std::cerr << "Warning: Atom number mismatch!" << std::endl;
      }

      if (j < pdbs.size()-1){
        delete mol;
      }
    }
  }

  if (trjout.is_open()){
    //Re-write header in case of any changes
    ftrjout->writeHeader(trjout);
    trjout.close();
    if (ftrjout != NULL){
      delete ftrjout;
    }
  }

  if (mol != NULL){
    delete mol;
  }

  return 0;
}
