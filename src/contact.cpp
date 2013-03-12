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
  std::cerr << "Usage:   contact [-options] <-pdb PDBFILE> <TRAJfile(s)>" << std::endl;
  std::cerr << "Options: [-sel selection]" << std::endl;
	std::cerr << "         [-add distance] [-scale value]" << std::endl;
	std::cerr << "         [-contact file]" << std::endl;
  std::cerr << "         [-skip frames] [-start frame]" << std::endl;
  exit(0);
}

int main (int argc, char **argv){

  int i;
	unsigned int j;
  std::vector<std::string> trajs;
  Molecule *mol;
	Molecule *cmol;
  std::string pdb;
  std::string sel=""; //To match NATOM in trajectory
	std::string fcontact;
	double add=0;
	double scale=1.0;
  std::string currArg;
  std::ifstream trjin;
	std::ifstream conFile;
	std::istream* coninp;
	Trajectory *ftrjin;
  int skip=0;
  int start=0;
	std::vector<Molecule *> p1;
	std::vector<Molecule *> p2;
	std::vector<double> d;
	std::vector<double> dist;
	std::vector<bool> conFlag;
	double minD;
	double conCount;
	std::string line;
	std::vector<std::string> s;

  pdb.clear();
	fcontact.clear();
  mol=NULL;
	cmol=NULL;
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
      sel=currArg;
    }
		else if (currArg == "-add"){
			currArg=argv[++i];
			std::stringstream(currArg) >> add;
		}
		else if (currArg == "-scale"){
			currArg=argv[++i];
			std::stringstream(currArg) >> scale;
		}
		else if (currArg == "-contact"){
			currArg=argv[++i];
			fcontact=currArg;
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
    if (sel.length() > 0){
      mol->select(sel); //Selection should match trajectory
      mol=mol->clone(true, false); //Delete original molecule after cloning
    }
  }

	//Set Up Contacts
	if (fcontact.length() > 0){
		conFile.open(fcontact.c_str(), std::ios::in);
		coninp=&conFile;
		while (coninp->good() && !(coninp->eof())){
			getline(*coninp, line);
			s.clear();
			s=Misc::split(line, " ", false);
			if (s.size() == 3){
				mol->select(s.at(0));
				cmol=mol->copy();
				p1.push_back(cmol);
				cmol=p1.at(0);
				mol->select(s.at(1));
        cmol=mol->copy();
	      p2.push_back(cmol);
				std::stringstream(s.at(2)) >> minD;
				d.push_back((minD+add)*scale);
			}
		}
		mol->selAll(); //Select all
		if (p1.size() != p2.size() || p1.size() != d.size()){
			std::cerr << "Error: Contact file input mismatch" << std::endl;
		}
	}
	else{
		
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
					//Get desired contacts
					dist.clear();
					conFlag.clear();
					conCount=0;
					for(j=0; j< p1.size(); j++){
						dist.push_back(Analyze::distance(p1.at(j), p2.at(j)));
						if (dist.at(j) <= d.at(j)){
							conFlag.push_back(true);
							conCount++;
						}
						else{
							conFlag.push_back(false);
						}
					}
					std::cout << conCount << " / " << conFlag.size();
					for (j=0; j< p1.size(); j++){
						std::cout << " " << conFlag.at(j) << " " << dist.at(j);
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
