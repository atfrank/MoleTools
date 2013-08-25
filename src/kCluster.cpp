//Sean M. Law

#include "Molecule.hpp"

#include <iostream>
#include <fstream>
#include <cstdlib>

void usage(){
  exit(0);
}

int main (int argc, char **argv){


  int i;
  std::string pdb;
  std::string currArg;
  std::string sel;

  pdb.clear();

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
      pdb=currArg;
    }
  }

  if (pdb.length() == 0){
    std::cerr << std::endl << "Error: Please provide an input file" << std::endl << std::endl;
    usage();
  }

  Molecule *mol=Molecule::readPDB(pdb);

  if (sel.length() >0){
    mol->select(sel);
  }

  Molecule *cmol=mol->clone();
  

  return 0;
}
