//Sean M. Law
//Aaron T. Frank
    
/*
This file is part of MoleTools.

MoleTools is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MoleTools is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MoleTools.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "Molecule.hpp"
#include "Residue.hpp"
#include "Atom.hpp"
#include "Coor.hpp"
#include "Misc.hpp"

#include <iostream>
#include <cstdlib>
#include <fstream>



void usage (){
  std::cerr << std::endl;
  std::cerr << "Usage:   modPDB [-options] <PDBfile>" << std::endl;
  std::cerr << "Options: [-model num]" << std::endl;
  std::cerr << "         [-renumber offset]" << std::endl;
  std::cerr << "         [-match matchPDB]" << std::endl;
  std::cerr << "         [-sel selection]" << std::endl;  
  std::cerr << "         [-fit fitPDB] [-fitsel selection]" << std::endl;
  std::cerr << "         [-data dataFile] [-bfac|occup] [-residue|atom] [-reset]" << std::endl;
  std::cerr << "         [-translate dx dy dz]" << std::endl;
  std::cerr << "         [-outsel selection]" << std::endl;
  std::cerr << "         [-rotate r1c1 r1c2 r1c3 r2c1 r2c2 r2c3 r3c1 r3c2 r3c3]" << std::endl;
  std::cerr << "         [-center] [-censel selection]" << std::endl;
  std::cerr << "         [-format type] [-chains]" << std::endl;
  std::cerr << "         [-nohetero]" << std::endl;
  std::cerr << "         [-cat PDBfile] [-catsel selection]" << std::endl;
  std::cerr << "         [-his]" << std::endl;
  std::cerr << std::endl << std::endl;
  exit(0);
}

int main (int argc, char **argv){

  int i;
  unsigned int j;
  int model=0;
  int residue_offset=0;
  bool renumber = false;
  bool match = false;
  std::string pdb;
  std::string match_pdb;
  std::string currArg;
  std::string sel;
  std::string fitsel=""; //Default = Fit All
  bool fit=false;
  std::string censel="";
  bool center=false;
  double dx, dy, dz;
  Coor dxyz;
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
  std::string datafile="";
  bool chnFlag;
  bool hetFlag;
  std::vector<std::string> catpdbs;
  bool hisFlag;
  bool data=false;
  bool data_bfac=true;
  bool data_occup=false;
  bool data_residue=true;
  bool data_atom=false;
  bool data_reset=false;

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
  hisFlag=false;

  for (i=1; i<argc; i++){
    currArg=argv[i];
    if (currArg.compare("-h") == 0 || currArg.compare("-help") == 0){
      usage();
    }
    else if (currArg.compare("-model") == 0){
      currArg=argv[++i];
      std::stringstream(currArg) >> model; //atoi
    }
    else if (currArg.compare("-renumber") == 0){
      currArg=argv[++i];
      renumber=true;
      std::stringstream(currArg) >> residue_offset; //atoi
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
    else if (currArg.compare("-match") == 0){
      currArg=argv[++i];
      match_pdb=currArg;
      match=true;
    }
    else if (currArg.compare("-data") == 0){
      currArg=argv[++i];
      datafile=currArg;
      data=true;
    }
    else if (currArg.compare("-bfac") == 0){
      data_bfac=true;
    }
    else if (currArg.compare("-occup") == 0){
      data_occup=true;
      data_bfac=false;
    }    
    else if (currArg.compare("-residue") == 0){
      data_residue=true;      
    }    
    else if (currArg.compare("-atom") == 0){
      data_atom=true;
      data_residue=false;
    }    
    else if (currArg.compare("-reset") == 0){
      data_reset=true;
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
      dxyz=Coor(dx, dy, dz);
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
      format=Misc::toupper(currArg);
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
    else if (currArg.compare("-his") == 0){
      hisFlag=true;
    }
    else{
      pdb=currArg;
    }
  }

  if (pdb.length() == 0){
    std::cerr << std::endl << "Error: Please provide an input PDB file" << std::endl << std::endl;
    usage();
  }

//  std::clock_t start;
//  double duration;
//  start=std::clock();

  mol=Molecule::readPDB(pdb, model, format, hetFlag);

  if (chnFlag == true){
    mol->addMissingChainIds();
  }
  if (hisFlag == true){
    mol->renameHis();
  }
  
  if (sel.length() >0){
    mol->select(sel);
    mol=mol->clone(true, false); //Clone and delete original
  }

//  duration = (std::clock() - start) / (double) CLOCKS_PER_SEC;
//  std::cerr << duration << std::endl;

  if (catpdbs.size() > 0){
    for (j=0; j< catpdbs.size(); j++){
      catmol=Molecule::readPDB(catpdbs.at(j), model, format, hetFlag);
      if (hisFlag == true){
        catmol->renameHis();
      }
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
      if (hisFlag == true){
        fitmol->renameHis();
      }
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

  else if (renumber){
  	Atom * atm;
  	int atom_num = 1;
  	mol->selAll();
  	for (j=0; j< mol->getAtmVecSize(); j++){
  		atm=mol->getAtom(j);
  		atm->setResId(atm->getResId()+residue_offset);
  		atm->setAtmNum(atom_num);
  		atom_num ++;
  	}
  }
  else if (match){
  	Molecule *matchmol;
  	matchmol=Molecule::readPDB(match_pdb, model, format, hetFlag);
  	outsel = ":.CA+CB+C+HA+H+N";
  	if ( matchmol->getResVecSize() == mol->getResVecSize()){
  		for (j=0; j< mol->getResVecSize(); j++){ 
  			mol->renameRes(mol->getResidue(j)->getAtom(0)->getResId(),mol->getResidue(j)->getResName(),matchmol->getResidue(j)->getResName());
  		}
  	}
  	else {
  		std::cerr << "Warning: match_pdb : " << match_pdb << " and pdb: " << pdb << " do not contain same number of residues" << std::endl;
  	}
  }  
  else if(data){
    std::string line;
    std::vector<std::string> s;  
		if (datafile.length() > 0){
		  std::ifstream datainp(datafile);
			while ( std::getline(datainp, line) ){
				Misc::splitStr(line, " ", s, true);
				if (s.size() >= 2){
					for (unsigned int j=0; j< mol->getAtmVecSize(); j++){
						if (data_residue){
							if (mol->getAtom(j)->getResId()==atoi(Misc::trim(s.at(0)).c_str())){
								if (data_occup){
									mol->getAtom(j)->setOccu(atof(Misc::trim(s.at(1)).c_str()));
								} 
								else {
									mol->getAtom(j)->setBFac(atof(Misc::trim(s.at(1)).c_str()));							
								}
							} else if(data_reset){
								if (data_occup){
									mol->getAtom(j)->setOccu(0.0);
								} else {
									mol->getAtom(j)->setBFac(0.0);
								}
							}
						}
						if (data_atom){
							if (mol->getAtom(j)->getAtmNum()==atoi(Misc::trim(s.at(0)).c_str())){
								if (data_occup){
									mol->getAtom(j)->setOccu(atof(Misc::trim(s.at(1)).c_str()));
								} 
								else {
									mol->getAtom(j)->setBFac(atof(Misc::trim(s.at(1)).c_str()));							
								}
							} else if(data_reset){
								if (data_occup){
									mol->getAtom(j)->setOccu(0.0);
								} else {
									mol->getAtom(j)->setBFac(0.0);
								}
							}
						}
					}						
				}					
			}
		} 	
  }
  else{
    //Do nothing
  }

  if (outsel.length() >0){
    mol->select(outsel);
  }

  mol->writePDB();

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
