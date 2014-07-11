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

#include "Cluster.hpp"

#include <iostream>

void usage(){
  std::cerr << std::endl << std::endl;
  std::cerr << "Usage:   cluster [-options] <DISTfile>" << std::endl;
  std::cerr << "Options: [-dp] [-k]" << std::endl;
  exit(0);
}

int main (int argc, char **argv){

  int i;
  std::string finp;
  std::string currArg;
  Cluster* clstr;

  finp.clear();
  clstr=NULL;

  for (i=1; i<argc; i++){
    currArg=argv[i];
    if (currArg.compare("-h") == 0 || currArg.compare("-help") == 0){
      usage();
    }
    else if (currArg.compare("-dp") == 0){
      clstr=new DPCluster;
    }
    else if (currArg.compare(0,1,"-") == 0){
      std::cerr << "Warning: Skipping unknown option \"" << currArg << "\"" << std::endl;
    }
    else{
      finp=currArg;
    }
  }

  if (clstr == NULL){
    clstr=new DPCluster;
  }

  clstr->readDMatrix(finp);

  return 0;
}
