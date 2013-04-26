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
  std::cerr << "Options: [-sel selection] [-cog selection]" << std::endl;
  std::cerr << "         [-dsel selection selection ] [-dist selection selection]" << std::endl;
	std::cerr << "         [-tsel selection selection selection] [-angle selection selection selection]" << std::endl;
	std::cerr << "         [-qsel selection selection selection selection]" << std::endl;
	std::cerr << "         [-dihedral selection selection selection selection]" << std::endl;
	std::cerr << "         [-fit refPDB] [-fitsel selection]" << std::endl;
	std::cerr << "         [-rmsd selection]" << std::endl;
  std::cerr << "         [-skip frames] [-start frame]" << std::endl;
  exit(0);
}

int main (int argc, char **argv){

  int i;
	unsigned int j;
  std::vector<std::string> trajs;
  Molecule *mol;
	Molecule *fitmol; //For fitting only
	Molecule *refmol; //For RMSD calculation only
  std::string pdb;
	std::string refpdb;
	bool fit=false;
	std::string fitsel=":.";
	std::string fcontact;
  std::string currArg;
  std::ifstream trjin;
	Trajectory *ftrjin;
  int skip=0;
  int start=0;
	std::string line;
	std::vector<std::string> s;
	Vector xyz;

	std::vector<Analyze *> analyses;
	Analyze *anin;

  pdb.clear();
	fcontact.clear();
  mol=NULL;
	fitmol=NULL;
	refmol=NULL;
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
    else if (currArg == "-sel" || currArg == "-nsel" || currArg == "-cog"){
			anin=new Analyze;
      anin->setType("quick");
      currArg=argv[++i];
			anin->addSel(currArg);
			analyses.push_back(anin);
    }
    else if (currArg == "-dsel" || currArg == "-dist" || currArg == "-distance"){
			anin=new Analyze;
      anin->setType("quick");
      currArg=argv[++i];
			anin->addSel(currArg);
      currArg=argv[++i];
			anin->addSel(currArg);
			analyses.push_back(anin);
    }
		else if (currArg == "-tsel" || currArg == "-angle"){
			anin=new Analyze;
      anin->setType("quick");
			currArg=argv[++i];
			anin->addSel(currArg);
			currArg=argv[++i];
			anin->addSel(currArg);
			currArg=argv[++i];
			anin->addSel(currArg);
			analyses.push_back(anin);
		}
		else if (currArg == "-qsel" || currArg == "-dihedral"){
			anin=new Analyze;
      anin->setType("quick");
      currArg=argv[++i];
			anin->addSel(currArg);
      currArg=argv[++i];
			anin->addSel(currArg);
      currArg=argv[++i];
			anin->addSel(currArg);
			currArg=argv[++i];
			anin->addSel(currArg);
			analyses.push_back(anin);
		}
		else if (currArg == "-fit"){
			fit=true;
			currArg=argv[++i];
			refpdb=currArg;
		}
		else if (currArg == "-fitsel"){
			fit=true;
			currArg=argv[++i];
			fitsel=currArg;
		}
		else if (currArg == "-rmsd"){
			anin=new Analyze;
			anin->setType("rmsd");
			currArg=argv[++i];
			anin->addSel(currArg);
      analyses.push_back(anin);
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

  if (fit == true){
    if (refpdb.length() > 0){
      refmol=Molecule::readPDB(refpdb);
      refmol->select(fitsel);
      refmol=refmol->clone(true, false); //Clone and delete original
    }
		//Else, use input PDB as reference, see below
  }

	if (pdb.length() == 0){
		std::cerr << std::endl << "Error: Please provide a PDB file via \"-pdb\"" << std::endl;
    usage();
	}
  else {
    mol=Molecule::readPDB(pdb);
		for (j=0; j< analyses.size(); j++){
			analyses.at(j)->setupMolSel(mol); //Make copies of mol from selection
			if (analyses.at(j)->getType() == "rmsd"){
				if (refmol == NULL){
					refmol=mol->clone();
				}
				analyses.at(j)->setupMolSel(refmol);
			}
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
					if (fit == true){
						ftrjin->getMolecule()->select(fitsel);
						ftrjin->getMolecule()->lsqfit(refmol);
						ftrjin->getMolecule()->selAll();
					}
					std::cout << ftrjin->getNPriv()*ftrjin->getTStepPS()/ftrjin->getNSavc()+i*ftrjin->getTStepPS();
					for (j=0; j< analyses.size(); j++){
						analyses.at(j)->runAnalysis();
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

	//delete mol;

  return 0;
}
