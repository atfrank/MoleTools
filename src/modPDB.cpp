//Sean M. Law

#include "Molecule.hpp"

#include <iostream>
#include <fstream>
#include <sstream> //string stream
#include <string>
#include <cstdlib>

#define MAXINPFILES 4096
#define MAXLINESIZE 4096

void usage (){
  std::cerr << std::endl;
  std::cerr << "Usage:   manipPDB [options] <PDBfile>" << std::endl;
  std::cerr << "Options: [-model num]" << std::endl;
  std::cerr << "         [-sel selection]" << std::endl;
  std::cerr << std::endl << std::endl;
  exit(0);
}

int main (int argc, char **argv){

  int i;
  int model=0;
  std::string pdb;
  std::string currArg;
  std::string sel;

  pdb.clear();
  
  for (i=1; i<argc; i++){
    currArg=argv[i];
    if (currArg == "-h" || currArg == "-help"){
      usage();
    }
    else if (currArg == "-model"){
      currArg=argv[++i];
      std::stringstream(currArg) >> model; //atoi
    }
    else if (currArg == "-sel" || currArg == "-nsel"){
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

  Molecule *mol=Molecule::readPDB(pdb, model);

  Molecule *cmol=mol->clone();
	
  cmol->deselAll();
  cmol->selAll();

  if (sel.length() >0){
    cmol->select(sel);
  }

//  mol->writePDB();
  cmol->writePDB();

	delete cmol;
	delete mol;

  return 0;
}
