//Sean M. Law

#include "Molecule.hpp"

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cstdlib>

void usage(){
  std::cerr << "Usage:   pcasso [-options] <pdbFile>" << std::endl;
  std::cerr << "Options: [-dssp dsspFile]" << std::endl;
	std::cerr << "         [-predict]" << std::endl;
//	std::cerr << "         [-trial]" << std::endl;
  std::cerr << std::endl;
  exit(0);
}

int main (int argc, char **argv){


  int i;
	unsigned int j;
  std::vector<std::string> pdbs;
  std::string currArg;
  std::string dssp;
	PcassoOutEnum out;

  dssp.clear();
	out=FEATURES;

  for (i=1; i<argc; i++){
    currArg=argv[i];
    if (currArg.compare("-h") == 0 || currArg.compare("-help") == 0){
      usage();
    }
    else if (currArg.compare("-dssp") == 0){
      currArg=argv[++i];
      dssp=currArg;
    }
		else if (currArg.compare("-predict") == 0 || currArg.compare("-prediction") == 0){
			out=PREDICT;
		}
    else{
      pdbs.push_back(currArg);
    }
  }

  if (pdbs.size() == 0){
    std::cerr << std::endl << "Error: Please provide an input file" << std::endl << std::endl;
    usage();
  }

  if (dssp.length() > 0){
    if (pdbs.size() > 1){
      std::cerr << "Warning: \"-dssp\" option can only be used with a single PDB input" << std::endl;
      std::cerr << "Warning: Only the first PDB structure was processed" << std::endl;
    }
    if (pdbs.size() > 0){
      //Only process first structure!
      Molecule *mol=Molecule::readPDB(pdbs.at(0));
      std::cerr << "Processing file \"" << pdbs.at(0) << "..." << std::endl;
			mol->pcasso(dssp, out); //Makes temporary clone with C-alpha only, and analyzes it
      delete mol;
    }
  }
  else{
	  for (j=0; j< pdbs.size(); j++){
  	  Molecule *mol=Molecule::readPDB(pdbs.at(j));
		  std::cerr << "Processing file \"" << pdbs.at(j) << "..." << std::endl;
			mol->pcasso("", out); //Makes temporary clone with C-alpha only, and analyzes it
		  delete mol;
    }
	}

  return 0;
}
