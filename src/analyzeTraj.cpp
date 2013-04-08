//Sean M. Law

#include "Trajectory.hpp"
#include "Molecule.hpp"
#include "Misc.hpp"
#include "Analyze.hpp"

#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <vector>

#define MAXINPFILES 4096
#define MAXLINESIZE 4096

void usage(){
  std::cerr << std::endl << std::endl;
  std::cerr << "Usage:   analyzeTraj [-options] <-pdb PDBFILE> <TRAJfile(s)>" << std::endl;
  std::cerr << "Options: [-sel selection]" << std::endl;
  std::cerr << "         [-dsel selection selection]" << std::endl;
	std::cerr << "         [-tsel selection selection selection]" << std::endl;
	std::cerr << "         [-qsel selection selection selection selection]" << std::endl;
  std::cerr << "         [-skip frames] [-start frame]" << std::endl;
  exit(0);
}

int main (int argc, char **argv){

  int i;
	unsigned int j;
  std::vector<std::string> trajs;
  Molecule *mol;
  Molecule *mol1;
  Molecule *mol2;
  Molecule *mol3;
  Molecule *mol4;
  std::string pdb;
  std::string sel1=""; //To match NATOM in trajectory
  std::string sel2="";
  std::string sel3="";
  std::string sel4="";
  int nsel=0; //Number of selections
	std::string fcontact;
  std::string currArg;
  std::ifstream trjin;
	Trajectory *ftrjin;
  int skip=0;
  int start=0;
	std::vector<Molecule *> p1;
	std::vector<Molecule *> p2;
	std::string line;
	std::vector<std::string> s;
	std::string analysis="";
	Vector xyz;

  pdb.clear();
	fcontact.clear();
  mol=NULL;
  mol1=NULL;
  mol2=NULL;
  mol3=NULL;
  mol4=NULL;
	ftrjin=NULL;

  for (i=1; i<argc; i++){
    currArg=argv[i];
    if (currArg == "-h" || currArg == "-help"){
      usage();
    }
    else if (currArg == "-pdb"){
      currArg=argv[++i];
      pdb=currArg;
    }
    else if (currArg == "-sel" || currArg == "-nsel"){
      currArg=argv[++i];
      sel1=currArg;
      nsel=1;
    }
    else if (currArg == "-dsel"){
      currArg=argv[++i];
      sel1=currArg;
      currArg=argv[++i];
      sel2=currArg;
      nsel=2;
    }
    else if (currArg == "-skip"){
      currArg=argv[++i];
      std::stringstream(currArg) >> skip;
    }
    else if (currArg == "-start"){
      currArg=argv[++i];
      std::stringstream(currArg) >> start;
      start--;
    }
    else{
      trajs.push_back(currArg);
    }
  }

  if (trajs.size() == 0){
    std::cerr << std::endl << "Error: Please provide an input trajectory file" << std::endl << std::endl;
    usage();
  }

	if (pdb.length() == 0){
		std::cerr << std::endl << "Error: Please provide a PDB file via \"-pdb\"" << std::endl;
    usage();
	}
  else {
    mol=Molecule::readPDB(pdb);
    if (sel1.length() > 0){
      mol->select(sel1); 
      mol1=mol->copy(); 
    }
    if (sel2.length() > 0){
      mol->select(sel2); 
      mol2=mol->copy(); 
    }
    if (sel3.length() > 0){
      mol->select(sel3); 
      mol3=mol->copy(); 
    }
    if (sel4.length() > 0){
      mol->select(sel4);
      mol4=mol->copy(); 
    } 
    mol->selAll();
  }

	//Process Trajectories

	for (j=0; j< trajs.size(); j++){
		trjin.open(trajs.at(j).c_str(), std::ios::binary);

		if (trjin.is_open()){
			ftrjin=new Trajectory;
      ftrjin->setMolecule(mol);

      if (ftrjin->findFormat(trjin) == true){
				ftrjin->readHeader(trjin);
        //Loop through desired frames
				for (i=start; i< ftrjin->getNFrame(); i=i+1+skip){
					ftrjin->readFrame(trjin, i);
					std::cout << i*ftrjin->getTStepPS();
					if (analysis == "SOME ANALYSIS"){

					}
					else{
						switch (nsel){
							case 1:
								xyz=Analyze::centerOfGeometry(mol1);
								std::cout << std::fixed;
								std::cout << std::setw(9) << std::right << std::setprecision(3) << xyz.x();
								std::cout << std::setw(9) << std::right << std::setprecision(3) << xyz.y();
								std::cout << std::setw(9) << std::right << std::setprecision(3) << xyz.z();
								break;
							case 2:
								std::cout << std::fixed;
								std::cout << std::setw(9) << std::right << std::setprecision(3) << Analyze::distance(Analyze::centerOfGeometry(mol1), Analyze::centerOfGeometry(mol2));
								break;
							case 3:
								std::cout << std::setw(9) << std::right << std::setprecision(3) << Analyze::angle(Analyze::centerOfGeometry(mol1), Analyze::centerOfGeometry(mol2), Analyze::centerOfGeometry(mol3));
								break;
							case 4:
								std::cout << std::setw(9) << std::right << std::setprecision(3) << Analyze::dihedral(Analyze::centerOfGeometry(mol1), Analyze::centerOfGeometry(mol2), Analyze::centerOfGeometry(mol3), Analyze::centerOfGeometry(mol4));
								break;
							default:
								break;
						}
					}
					std::cout << std::endl;
				}
			}
			else{
				std::cerr << "Warning: Skipping unknown trajectory format \"";
				std::cerr << trajs.at(j) << "\"" << std::endl;
			}
			if (ftrjin != NULL){
				delete ftrjin;
			}
		}
	
		trjin.close();
	}

	for (j=0; j< p1.size(); j++){
		delete p1.at(j);
	}
	for (j=0; j< p2.size(); j++){
		delete p2.at(j);
  }
	//delete mol;

  return 0;
}
