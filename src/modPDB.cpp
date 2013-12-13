//Sean M. Law

#include "Molecule.hpp"
#include "Vector.hpp"

#include <iostream>
#include <fstream>
#include <sstream> //string stream
#include <string>
#include <cstdlib>

void usage (){
  std::cerr << std::endl;
  std::cerr << "Usage:   modPDB [-options] <PDBfile>" << std::endl;
  std::cerr << "Options: [-model num]" << std::endl;
  std::cerr << "         [-sel selection]" << std::endl;
	std::cerr << "         [-fit fitPDB] [-fitsel selection]" << std::endl;
  std::cerr << "         [-translate dx dy dz]" << std::endl;
  std::cerr << "         [-outsel selection]" << std::endl;
	std::cerr << "         [-rotate r1c1 r1c2 r1c3 r2c1 r2c2 r2c3 r3c1 r3c2 r3c3]" << std::endl;
  std::cerr << "         [-center] [-censel selection]" << std::endl;
	std::cerr << "         [-format type] [-chains]" << std::endl;
  std::cerr << "         [-nohetero]" << std::endl;
	std::cerr << "         [-cat PDBfile] [-catsel selection]" << std::endl;
  std::cerr << std::endl << std::endl;
  exit(0);
}

int main (int argc, char **argv){

  int i;
	unsigned int j;
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
  std::string outsel="";
	Molecule *mol;
	Molecule *fitmol;
	Molecule *catmol;
	std::string catsel="";
	std::string format;
	bool chnFlag;
  bool hetFlag;
	std::vector<std::string> catpdbs;

  pdb.clear();
	mol=NULL;
	fitmol=NULL; //Stationary molecule
	catmol=NULL;
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
	format.clear();
	chnFlag=false;
  hetFlag=true; //Keep hetero atoms
  catpdbs.clear();

  for (i=1; i<argc; i++){
    currArg=argv[i];
    if (currArg.compare("-h") == 0 || currArg.compare("-help") == 0){
      usage();
    }
    else if (currArg.compare("-model") == 0){
      currArg=argv[++i];
      std::stringstream(currArg) >> model; //atoi
    }
    else if (currArg.compare("-sel") == 0 || currArg.compare("-nsel") == 0){
      currArg=argv[++i];
      sel=currArg;
    }
    else if (currArg.compare("-outsel") == 0){
      currArg=argv[++i];
      outsel=currArg;
    }
		else if (currArg.compare("-fitsel") == 0){
			currArg=argv[++i];
			fitsel=currArg;
			fit=true;
		}
		else if (currArg.compare("-fit") == 0){
			currArg=argv[++i];
			fitpdb=currArg;
			fit=true;
		}
    else if (currArg.compare("-center") == 0){
      center=true;
    }
    else if (currArg.compare("-censel") == 0){
      currArg=argv[++i];
      censel=currArg;
      center=true;
    }
    else if (currArg.compare("-translate") == 0){
      currArg=argv[++i];
      std::stringstream(currArg) >> dx;
      currArg=argv[++i];
      std::stringstream(currArg) >> dy;
      currArg=argv[++i];
      std::stringstream(currArg) >> dz;
      dxyz=Vector(dx, dy, dz);
      translate=true;
    }
		else if (currArg.compare("-rotate") == 0){
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
		else if (currArg.compare("-format") == 0){
			currArg=argv[++i];
			Misc::toupper(currArg);
			format=currArg;
		}
		else if (currArg.compare("-chains") == 0){
			chnFlag=true;
		}
    else if (currArg.compare("-nohetero") == 0){
      hetFlag=false;
    }
		else if (currArg.compare("-cat") == 0){
			currArg=argv[++i];
			catpdbs.push_back(currArg);
		}
		else if (currArg.compare("-catsel") == 0){
			currArg=argv[++i];
			catsel=currArg;
		}
    else{
      pdb=currArg;
    }
  }

  if (pdb.length() == 0){
    std::cerr << std::endl << "Error: Please provide an input PDB file" << std::endl << std::endl;
    usage();
  }

  mol=Molecule::readPDB(pdb, model, format, hetFlag);
	
  if (sel.length() >0){
    mol->select(sel);
		mol=mol->clone(true, false); //Clone and delete original
  }

	if (catpdbs.size() > 0){
		for (j=0; j< catpdbs.size(); j++){
			catmol=Molecule::readPDB(catpdbs.at(j), model, format, hetFlag);
			if (catsel.length() > 0){
				catmol->select(catsel);
				catmol=catmol->clone(true, false); //Clone and delete original
			}
			mol->cat(catmol,true,false); //Cat selected and delete original
			catmol=NULL;
		}
	}

	if (fit == true){
		if (fitpdb.length() > 0){
			fitmol=Molecule::readPDB(fitpdb, model, format, hetFlag);
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

  if (outsel.length() >0){
    mol->select(outsel);
  }
  mol->writePDB(chnFlag);

	if (mol != NULL){
		delete mol;
	}

	if (fitmol != NULL){
		delete fitmol;
	}

	if (catmol != NULL){
		delete catmol;
	}

  return 0;
}
