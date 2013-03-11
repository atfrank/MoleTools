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
	std::cerr << "         [-fit refPDB] [-fitsel selection]" << std::endl;
  std::cerr << std::endl << std::endl;
  exit(0);
}

int main (int argc, char **argv){

  int i;
  int model=0;
  std::string pdb;
  std::string currArg;
  std::string sel;
	std::string fitsel=":."; //Default = Fit All
	bool fit=false;
	std::string refpdb;
	Molecule *mol;
	Molecule *refmol;

  pdb.clear();
	mol=NULL;
	refmol=NULL; //Stationary molecule
  
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
		else if (currArg == "-fitsel"){
			currArg=argv[++i];
			fitsel=currArg;
			fit=true;
		}
		else if (currArg == "-fit"){
			currArg=argv[++i];
			refpdb=currArg;
			fit=true;
		}
    else{
      pdb=currArg;
    }
  }

  if (pdb.length() == 0){
    std::cerr << std::endl << "Error: Please provide an input PDB file" << std::endl << std::endl;
    usage();
  }

  mol=Molecule::readPDB(pdb, model);
  if (sel.length() >0){
    mol->select(sel);
		mol=mol->clone(true, false); //Clone and delete original
  }

	if (fit == true){
		if (refpdb.length() > 0){
			refmol=Molecule::readPDB(refpdb, model);
			refmol->select(fitsel);
			mol->select(fitsel);
	  	mol->lsqfit(refmol);
			mol->rmsd(refmol);
		}
		else {
			std::cerr << std::endl << "Error: Please provide a reference PDB file for fitting" << std::endl;
			usage();
		}
	}

  mol->writePDB();

	if (mol != NULL){
		delete mol;
	}

  return 0;
}
