//Sean M. Law

#include "Molecule.hpp"

#include <iostream>
#include <fstream>
#include <cstdlib>

void usage(){
  std::cerr << std::endl;
  std::cerr << "Usage:   mol22PDB [-options] <MOL2file>" << std::endl;
  std::cerr << "Options: [-sel selection]" << std::endl;
	std::cerr << "         [-format type] [-chains]" << std::endl;
//	std::cerr << "         [-warnings]" << std::endl;
  std::cerr << std::endl;
  exit(0);
}

int main (int argc, char **argv){


  int i;
  std::string mol2;
  std::string currArg;
  std::string sel;
	std::string format;
  bool chnFlag;
//	bool warnings;

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
			Misc::toupper(currArg);
			format=currArg;
		}
    else if (currArg.compare("-chains") == 0){
      chnFlag=true;
    }
    else{
      mol2=currArg;
    }
  }

  if (mol2.length() == 0){
    std::cerr << std::endl << "Error: Please provide an input file" << std::endl << std::endl;
    usage();
  }

  Molecule *mol=Molecule::readMol2(mol2, format);

  if (sel.length() >0){
    mol->select(sel);
  }

 	mol->writePDB(chnFlag); 

  return 0;
}
