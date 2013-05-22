//Sean M. Law

#include "Molecule.hpp"

#include <iostream>
#include <fstream>
#include <cstdlib>

#define MAXINPFILES 4096
#define MAXLINESIZE 4096

void usage(){
  std::cerr << std::endl;
  std::cerr << "Usage:   mol22PDB [-options] <MOL2file>" << std::endl;
  std::cerr << "Options: [-sel selection]" << std::endl;
	std::cerr << "         [-format type] [-chains]" << std::endl;
  std::cerr << std::endl;
  exit(0);
}

int main (int argc, char **argv){


  int i;
  std::string mol2;
  std::string currArg;
  std::string sel;
	std::string format="UNFORMATTED";
  bool chnFlag;

  mol2.clear();
  chnFlag=false;

  for (i=1; i<argc; i++){
    currArg=argv[i];
    if (currArg == "-h" || currArg == "-help"){
      usage();
    }
    else if (currArg == "-sel" || currArg == "-nsel"){
      currArg=argv[++i];
      sel=currArg;
    }
		else if (currArg == "-format"){
			currArg=argv[++i];
			Misc::toupper(currArg);
			format=currArg;
		}
    else if (currArg == "-chains"){
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

  Molecule *mol=Molecule::readMol2(mol2);

  if (sel.length() >0){
    mol->select(sel);
  }

 	mol->writePDB(true, true, chnFlag, format); 

  return 0;
}
