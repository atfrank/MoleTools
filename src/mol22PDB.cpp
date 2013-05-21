//Sean M. Law

#include "Molecule.hpp"

#include <iostream>
#include <fstream>
#include <cstdlib>

#define MAXINPFILES 4096
#define MAXLINESIZE 4096

void usage(){
  exit(0);
}

int main (int argc, char **argv){


  int i;
  std::string mol2;
  std::string currArg;
  std::string sel;

  mol2.clear();

  for (i=1; i<argc; i++){
    currArg=argv[i];
    if (currArg == "-h" || currArg == "-help"){
      usage();
    }
    else if (currArg == "-nsel"){
      currArg=argv[++i];
      sel=currArg;
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

 	mol->writePDB(true, true, true); 

  return 0;
}
