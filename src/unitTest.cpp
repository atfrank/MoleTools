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

#include <cstdlib>

void usage(){
  std::cerr << std::endl << std::endl;
  std::cerr << "Usage:   unitTest [-options] <file>" << std::endl;
  std::cerr << "Options: [-format type]" << std::endl;
  std::cerr << "         [-top file] [-prm file]" << std::endl;
  exit(0);
}

int main (int argc, char **argv){


  int i;
  std::string currArg;
	std::vector<std::string> pdbs;

	pdbs.clear();

  for (i=1; i<argc; i++){
    currArg=argv[i];
    if (currArg.compare("-h") == 0 || currArg.compare("-help") == 0){
      usage();
    }
		else if (currArg.compare(0,1,"-") == 0){
      std::cerr << "Warning: Skipping unknown option \"" << currArg << "\"" << std::endl;
    }
    else{
			pdbs.push_back(currArg);
    }
  }

	for (unsigned int j=0; j< pdbs.size(); j++){
		Molecule *mol=Molecule::readPDB(pdbs.at(j));

		mol->select(":.CA");
		Molecule *cmol=mol->copy(true);

		delete cmol;
		delete mol;
	}

  return 0;
}
