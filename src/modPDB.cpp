//Sean M. Law

#include "Molecule.hpp"
#include "Vector.hpp"

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
	std::cerr << "         [-fit fitPDB] [-fitsel selection]" << std::endl;
  std::cerr << "         [-translate dx dy dz]" << std::endl;
	std::cerr << "         [-rotate r1c1 r1c2 r1c3 r2c1 r2c2 r2c3 r3c1 r3c2 r3c3]" << std::endl;
  std::cerr << "         [-center] [-censel selection]" << std::endl;
  std::cerr << std::endl << std::endl;
  exit(0);
}

int main (int argc, char **argv){

  int i;
  int model=0;
  std::string pdb;
  std::string currArg;
  std::string sel;
	std::string fitsel=""; //Default = Fit All
	bool fit=false;
  std::string censel="";
  bool center=false;
  double dx, dy, dz;
  Vector dxyz;
  bool translate=false;
	double r1c1, r1c2, r1c3, r2c1, r2c2, r2c3, r3c1, r3c2, r3c3;
	bool rotate=false;
	std::string fitpdb;
	Molecule *mol;
	Molecule *fitmol;

  pdb.clear();
	mol=NULL;
	fitmol=NULL; //Stationary molecule
  dx=0.0;
  dy=0.0;
  dz=0.0;
	r1c1=0.0;
	r1c2=0.0;
	r1c3=0.0;
 	r2c1=0.0;
	r2c2=0.0;
	r2c3=0.0;
	r3c1=0.0;
	r3c2=0.0;
	r3c3=0.0;
  
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
			fitpdb=currArg;
			fit=true;
		}
    else if (currArg == "-center"){
      center=true;
    }
    else if (currArg == "-censel"){
      currArg=argv[++i];
      censel=currArg;
      center=true;
    }
    else if (currArg == "-translate"){
      currArg=argv[++i];
      std::stringstream(currArg) >> dx;
      currArg=argv[++i];
      std::stringstream(currArg) >> dy;
      currArg=argv[++i];
      std::stringstream(currArg) >> dz;
      dxyz=Vector(dx, dy, dz);
      translate=true;
    }
		else if (currArg == "-rotate"){
			currArg=argv[++i];
      std::stringstream(currArg) >> r1c1;
      currArg=argv[++i];
      std::stringstream(currArg) >> r1c2;
      currArg=argv[++i];
      std::stringstream(currArg) >> r1c3;
			currArg=argv[++i];
      std::stringstream(currArg) >> r2c1;
      currArg=argv[++i];
      std::stringstream(currArg) >> r2c2;
      currArg=argv[++i];
      std::stringstream(currArg) >> r2c3;
			currArg=argv[++i];
      std::stringstream(currArg) >> r3c1;
      currArg=argv[++i];
      std::stringstream(currArg) >> r3c2;
      currArg=argv[++i];
      std::stringstream(currArg) >> r3c3;
      rotate=true;
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
		if (fitpdb.length() > 0){
			fitmol=Molecule::readPDB(fitpdb, model);
      if (fitsel.length() > 0){
			  fitmol->select(fitsel);
			  mol->select(fitsel);
      }
      else{
        fitmol->selAll();
        mol->selAll();
      }
	  	mol->lsqfit(fitmol);
			//mol->rmsd(fitmol);
      mol->selAll();
      delete fitmol;
		}
		else {
			std::cerr << std::endl << "Error: Please provide a PDB file for fitting" << std::endl;
			usage();
		}
	}
  else if (center == true){
    if (censel.length() > 0){ 
      mol->select(censel);
    }
    else{
      mol->selAll();
    }
    mol->center();
    mol->selAll();
  }
  else if (translate == true){
    mol->translate(dxyz);
  }
	else if (rotate == true){
		mol->rotate(r1c1, r1c2, r1c3, r2c1, r2c2, r2c3, r3c1, r3c2, r3c3);
	}
  else{
    //Do nothing
  }

  mol->writePDB();

	if (mol != NULL){
		delete mol;
	}

  return 0;
}
