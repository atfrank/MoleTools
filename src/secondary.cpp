//Sean M. Law

#include "Molecule.hpp"
#include "SABA.hpp"

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cstdlib>

void usage(){
  exit(0);
}

int main (int argc, char **argv){


  int i;
	unsigned int j;
  std::vector<std::string> pdbs;
  std::string currArg;

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

  if (pdbs.size() == 0){
    std::cerr << std::endl << "Error: Please provide an input file" << std::endl << std::endl;
    usage();
  }

	for (j=0; j< pdbs.size(); j++){
  	Molecule *mol=Molecule::readPDB(pdbs.at(j));
		mol=SABA::getPseudoCenter(mol);
		delete mol;
	}

  return 0;
}
