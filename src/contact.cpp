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

void usage(){
  std::cerr << std::endl << std::endl;
  std::cerr << "Usage:   contact [-options] <-pdb PDBFILE> <TRAJfile(s)>" << std::endl;
  std::cerr << "Options: [-sel selection]" << std::endl;
	std::cerr << "         [-add distance] [-scale value]" << std::endl;
	std::cerr << "         [-ref file]" << std::endl;
  std::cerr << "         [-skip frames] [-start frame]" << std::endl;
  std::cerr << "         [-prob] [-silent]" << std::endl;
	std::cerr << "         [-verbose]" << std::endl;
  exit(0);
}

int main (int argc, char **argv){

  int i;
	unsigned int j, k;
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
  bool probFlag;
  std::vector<double> probs;
  unsigned int ndata;
  std::vector<std::string> sels;
	bool silentFlag;
	bool verboseFlag;


  pdb.clear();
	fcontact.clear();
  mol=NULL;
	cmol=NULL;
	ftrjin=NULL;
  probFlag=false;
  ndata=0;
	silentFlag=false;
	verboseFlag=false;

  for (i=1; i<argc; i++){
    currArg=argv[i];
    if (currArg.compare("-h") == 0 || currArg.compare("-help") == 0){
      usage();
    }
    else if (currArg.compare("-pdb") == 0){
      currArg=argv[++i];
      pdb=currArg;
    }
    else if (currArg.compare("-sel") == 0 || currArg.compare("-nsel") == 0){
      currArg=argv[++i];
      sel=currArg;
    }
		else if (currArg.compare("-add") == 0){
			currArg=argv[++i];
			std::stringstream(currArg) >> add;
		}
		else if (currArg.compare("-scale") == 0){
			currArg=argv[++i];
			std::stringstream(currArg) >> scale;
		}
		else if (currArg.compare("-ref") == 0){
			currArg=argv[++i];
			fcontact=currArg;
		}
    else if (currArg.compare("-skip") == 0){
      currArg=argv[++i];
      std::stringstream(currArg) >> skip;
    }
    else if (currArg.compare("-start") == 0){
      currArg=argv[++i];
      std::stringstream(currArg) >> start;
      start--;
    }
    else if (currArg.compare("-prob") == 0){
      probFlag=true;
    }
		else if (currArg.compare("-silent") == 0){
			silentFlag=true;		
		}
		else if (currArg.compare("-verbose") == 0){
			verboseFlag=true;
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
			Misc::splitStr(line, " \t", s, false);
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
        if (probFlag == true){
          sels.push_back(line);
        }
			}
		}
		mol->selAll(); //Select all
		if (p1.size() != p2.size() || p1.size() != d.size()){
			std::cerr << "Error: Contact file input mismatch" << std::endl;
		}
	}
	else{
		
	}

  if (probFlag == true){
    probs.resize(p1.size(),0);
  }

	//Process Trajectories
	for (j=0; j< trajs.size(); j++){
		trjin.open(trajs.at(j).c_str(), std::ios::binary);

		if (trjin.is_open()){
			ftrjin=new Trajectory;
      ftrjin->setMolecule(mol);
			
      if (ftrjin->findFormat(trjin) == true){
				if (verboseFlag == true){
					std::cerr << "Processing file \"" << trajs.at(j).c_str() << "\"..." << std::endl;
				}
				ftrjin->readHeader(trjin);
        //Loop through desired frames
				for (i=start; i< ftrjin->getNFrame(); i=i+1+skip){
					ftrjin->readFrame(trjin, i);
					//Get desired contacts
					dist.clear();
					conFlag.clear();
					conCount=0;
					for(k=0; k< p1.size(); k++){
						dist.push_back(Analyze::distance(p1.at(k), p2.at(k)));
						if (dist.at(k) <= d.at(k)){
							conFlag.push_back(true);
							conCount++;
              if (probFlag == true){
                probs.at(k)+=1.0;
              }
						}
						else{
							conFlag.push_back(false);
						}
					}
					if (silentFlag == false){
						std::cout << conCount << " / " << conFlag.size();
						for (k=0; k< p1.size(); k++){
							std::cout << " " << conFlag.at(k) << " " << dist.at(k);
						}
						std::cout << std::endl;
					}
          if (probFlag == true){
            ndata++;
          }
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

  if (probFlag == true){
    for (k=0; k< probs.size(); k++){
      std::cerr << sels.at(k) << " Probability = " << probs.at(k)/ndata << std::endl;
    }
  }

  return 0;
}
