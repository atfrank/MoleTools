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
#include "Misc.hpp"

void usage(){
  std::cerr << std::endl;
  std::cerr << "Usage:   mol22PDB [-options] <MOL2files>" << std::endl;
  std::cerr << "Options: [-sel selection]" << std::endl;
  std::cerr << "         [-format type] [-chains]" << std::endl;
//  std::cerr << "         [-warnings]" << std::endl;
  std::cerr << std::endl;
  exit(0);
}

int main (int argc, char **argv){


  int i;
  std::vector<std::string> mol2;
  std::string currArg;
  std::string sel;
  std::string format;
  bool chnFlag;
//  bool warnings;

  mol2.clear();
  chnFlag=false;
  format.clear();

  for (i=1; i<argc; i++){
    currArg=argv[i];
    if (currArg.compare("-h") == 0 || currArg.compare("-help") == 0){
      usage();
    }
    else if (currArg.compare("-sel") == 0 || currArg.compare("-nsel") == 0){
      currArg=argv[++i];
      sel=currArg;
    }
    else if (currArg.compare("-format") == 0){
      currArg=argv[++i];
      format=Misc::toupper(currArg);
    }
    else if (currArg.compare("-chains") == 0){
      chnFlag=true;
    }
    else if (currArg.compare(0,1,"-") == 0){
      std::cerr << "Warning: Skipping unknown option \"" << currArg << "\"" << std::endl;
    }
    else{
      mol2.push_back(currArg);
    }
  }

  if (mol2.size() == 0){
    std::cerr << std::endl << "Error: Please provide an input file" << std::endl << std::endl;
    usage();
  }

  Molecule *mol=new Molecule;

  for (unsigned int j=0; j< mol2.size(); j++){
    mol->cat(Molecule::readMol2(mol2.at(j), format));
  }

  if (sel.length() >0){
    mol->select(sel);
  }

  mol->writePDB(chnFlag);

  mol=mol->clone();
  for (unsigned int j=0; j< mol->getNAtom(); j++){
//    std::cerr << mol->getAtom(j)->getAtmType() << std::endl;
  }

  delete mol;

  return 0;
}
