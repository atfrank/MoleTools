//Sean M. Law

#include "Trajectory.hpp"
#include "Molecule.hpp"

#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>

#define MAXINPFILES 4096
#define MAXLINESIZE 4096

void usage(){
  std::cerr << std::endl << std::endl;
  std::cerr << "Usage:   pdb2Traj [-options] <PDBfile(s)>" << std::endl;
  std::cerr << "         [-out TRAJname] [-outsel selection]" << std::endl;
  exit(0);
}

int main (int argc, char **argv){

  int i;
	unsigned int j;
  std::string fout;
  Molecule *mol;
  std::vector<std::string> pdbs;
  std::string currArg;
  std::string outsel="";
  std::ofstream trjout;
	Trajectory *ftrjout;

  pdbs.clear();
  mol=NULL;
	ftrjout=NULL;

  for (i=1; i<argc; i++){
    currArg=argv[i];
    if (currArg == "-h" || currArg == "-help"){
      usage();
    }
    else if (currArg == "-out"){
      currArg=argv[++i];
      fout=currArg;
    }
    else if (currArg == "-outsel"){
      currArg=argv[++i];
      outsel=currArg;
    }
    else{
      pdbs.push_back(currArg);
    }
  }

  if (pdbs.size() == 0){
    std::cerr << std::endl << "Error: Please provide an input trajectory file" << std::endl << std::endl;
    usage();
  }
  else if (fout.length() == 0){
    std::cerr << std::endl << "Error: Please specify an output trajectory via \"-out\"" << std::endl; 
  }
	else{
		trjout.open(fout.c_str(), std::ios::binary);
    ftrjout=new Trajectory;
		ftrjout->setDefaultHeader();
		ftrjout->setNFrame(static_cast<int>(pdbs.size()));

		for (j=0; j< pdbs.size(); j++){
			mol=Molecule::readPDB(pdbs.at(j));
			if (outsel.length() > 0){
				mol->select(outsel);
				mol=mol->clone(true, false);
			}
			
			if (j==0){
				ftrjout->setNAtom(static_cast<int>(mol->getNAtom()));
				ftrjout->writeHeader(trjout);
			}

			if (ftrjout->getNAtom() == static_cast<int>(mol->getNAtom())){
   			ftrjout->setMolecule(mol);
	
				if (trjout.is_open()){
     			ftrjout->writeFrame(trjout);
   			}
			}
			else{
				ftrjout->setNFrame(ftrjout->getNFrame()-1);
				std::cerr << "Warning: Atom number mismatch!" << std::endl;
			}

			delete mol;
		}
	}

  if (trjout.is_open()){
    //Re-write header in case of any changes
    ftrjout->writeHeader(trjout);
    trjout.close();
		if (ftrjout != NULL){
			delete ftrjout;
		}
  }

  return 0;
}