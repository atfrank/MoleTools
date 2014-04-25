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
#include "Chain.hpp"
#include "Residue.hpp"
#include "Atom.hpp"
#include "Misc.hpp"

#include <iostream>
#include <fstream>

void usage(){
  exit(0);
}

int main (int argc, char **argv){


  int i;
  std::vector<std::string> pdbs;
  std::string currArg;
  std::string flist;
  std::fstream listFile;
  std::istream* listinp;
  std::string line;
  Chain* chn;
  Residue *res;
  Atom *atm;

  pdbs.clear();
  flist.clear();
  chn=NULL;
  res=NULL;
  atm=NULL;

  for (i=1; i<argc; i++){
    currArg=argv[i];
    if (currArg.compare("-h") == 0 || currArg.compare("-help") == 0){
      usage();
    }
    else if (currArg.compare("-list") == 0){
      currArg=argv[++i];
      flist=currArg;
    }
    else if (currArg.compare(0,1,"-") == 0){
      std::cerr << "Warning: Skipping unknown option \"" << currArg << "\"" << std::endl;
    }
    else{
      pdbs.push_back(currArg);
    }
  }

  if (pdbs.size() == 0 && flist.length() == 0){
    std::cerr << std::endl << "Error: Please provide an input file" << std::endl << std::endl;
    usage();
  }

  if (flist.length() > 0){
    listFile.open(flist.c_str(), std::ios::in);
    listinp=&listFile;
    while (listinp->good() && !(listinp->eof())){
      getline(*listinp, line);
      if (line.length() > 0){
        pdbs.push_back(line);
      }
    }
  }

  for (unsigned int j=0; j< pdbs.size(); j++){
    Molecule *mol=Molecule::readPDB(pdbs.at(j));
    std::cerr << "Processing file \"" << pdbs.at(j) << "\"..." << std::endl;

    mol->select(":protein+unk.backbone", false);
    Molecule *cmol=mol->copy(true);

    if (cmol->getNAtom() != 0){
      bool done=false;
      for (unsigned int c=0; c< cmol->getChnVecSize(); c++){
        chn=cmol->getChain(c);
        for (unsigned int r=0; r< chn->getResVecSize(); r++){
          if (r == 0 || r == chn->getResVecSize()){
            continue;
          }
          res=chn->getResidue(r);
          unsigned int count=0;
          for (unsigned int a=0; a< res->getAtmVecSize(); a++){
            atm=res->getAtom(a);
            std::string atmname=Misc::trim(atm->getAtmName());
            if (atmname.compare("CA") == 0){
              count++;
            }
            else if (atmname.compare("N") == 0){
              count++;
            }
            else if (atmname.compare("C") == 0){
              count++;
            }
            else if (atmname.compare("O") == 0){
              count++;
            }
            else if (atmname.compare(0,2,"OT") == 0){
              count++;
            }
            else{
              //Do Nothing
            }
            if (count >= 4){
              break;
            }
          }
          if (count < 4){
            done=true;
            break;
          }
        }
        if (done == true){
          break;
        }
      }

      std::cout << pdbs.at(j) << " " << done << " " << cmol->getYear() << " ";
      if (done == true && atm != NULL){
        std::cout << atm->getChainId() << ":";
        std::cout << atm->getResName() << atm->getResId() << ". ";
      }
      else{
        std::cout << ":. ";
      }
      if (mol->getExp().length() > 0){
        std::cout << Misc::replace(mol->getExp(), " ", "_", true) << " ";
      } 
      else{
        std::cout << "N/A ";
      }
      std::cout << std::endl;
    }
    else{
      std::cout << "# " << pdbs.at(j) << " " << cmol->getNAtom() << std::endl;
    }
    
    delete mol;
    delete cmol;
  }

  return 0;
}
