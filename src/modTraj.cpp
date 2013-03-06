//Sean M. Law

#include "Trajectory.hpp"
#include "Molecule.hpp"

#include <iostream>
#include <fstream>
#include <cstdlib>

#define MAXINPFILES 4096
#define MAXLINESIZE 4096

void usage(){
  std::cerr << std::endl << std::endl;
  std::cerr << "Usage:   modTraj [-options] <TRAJfile(s)>" << std::endl;
  std::cerr << "Options: [-pdb PDBfile] [-sel selection]" << std::endl;
  std::cerr << "         [-fit | -recenter] [-fitsel selection] [-recsel selection]" << std::endl;
  std::cerr << "         [-out TRAJname] [-outsel selection]" << std::endl;
  std::cerr << "         [-skip frames]" << std::endl;
	std::cerr << "         [-verbose] [-show]" << std::endl;
  exit(0);
}

int main (int argc, char **argv){


  int i;
	unsigned int j;
  std::vector<std::string> trajs;
  std::string out;
  Molecule *mol;
  Molecule *outmol;
  Molecule *fitmol;
  Molecule *recmol;
  std::string pdb;
  std::string sel=""; //To match NATOM in trajectory
  std::string currArg;
  bool fit=false;
  std::string fitsel=":.CA+P";
  bool recenter=false;
  std::string recsel=":.CA+P";
  std::string outsel="";
  std::ifstream trjin;
  std::ofstream trjout;
	Trajectory *ftrjin;
	Trajectory *ftrjout;
  bool show=false;

  pdb.clear();
  mol=NULL;
  outmol=NULL;
  fitmol=NULL;
  recmol=NULL;
	ftrjin=NULL;
	ftrjout=NULL;

  for (i=1; i<argc; i++){
    currArg=argv[i];
    if (currArg == "-h" || currArg == "-help"){
      usage();
    }
    else if (currArg == "-pdb"){
      currArg=argv[++i];
      pdb=currArg;
    }
    else if (currArg == "-sel"){
      currArg=argv[++i];
      sel=currArg;
    }
    else if (currArg == "-fit"){
      fit=true;
    }
    else if (currArg == "-fitsel"){
      currArg=argv[++i];
      fitsel=currArg;
      fit=true;
    }
    else if (currArg == "-recenter"){
      recenter=true;
    }
    else if (currArg == "recsel"){
      currArg=argv[++i];
      recsel=currArg;
      recenter=true;
    }
    else if (currArg == "-out"){
      currArg=argv[++i];
      out=currArg;
    }
    else if (currArg == "-outsel"){
      currArg=argv[++i];
      outsel=currArg;
    }
    else if (currArg == "-show"){
      show=true;
    }
    else{
      trajs.push_back(currArg);
    }
  }

  if (trajs.size() == 0){
    std::cerr << std::endl << "Error: Please provide an input trajectory file" << std::endl << std::endl;
    usage();
  }
  else if (outsel.length() > 0 && out.length() == 0){
    std::cerr << std::endl << "Error: Please specify an output trajectory via \"-out\" ";
    std::cerr << std::endl << "when using option \"-outsel\"" << std::endl;
  }

  if (fit == true && recenter == true){
    std::cerr << std::endl << "Error: Options \"-fit\" and \"-recenter\" cannot be used simultaneously";
    std::cerr << std::endl;
    usage();
  }

  if ((outsel.length() > 0 || fit == true || recenter == true) && pdb.length() == 0){
    std::cerr << std::endl << "Error: Please provide a PDB file via \"-pdb\" ";
    std::cerr << "in order to make a selection" << std::endl;
    usage();
  }
  else if (pdb.length() > 0){
    mol=Molecule::readPDB(pdb);
    if (sel.length() > 0){
      mol->select(sel); //Selection should match trajectory
      mol=mol->clone(true, false); //Delete original molecule after cloning
    }

    if (outsel.length() > 0){
      mol->select(outsel);
      outmol=mol->copy();
    }
    if (fit == true){
      mol->select(fitsel);
      fitmol=mol->copy();
    }
    if (recenter == true){
      mol->select(recsel);
      recmol=mol->copy();
    }
  }
  else{
    //Do Nothing
  }

  if (out.length() > 0){
    trjout.open(out.c_str(), std::ios::binary);
		ftrjout=new Trajectory;
    if (outmol != NULL){
      ftrjout->setMolecule(outmol);
    }
  }

	for (j=0; j< trajs.size(); j++){
		trjin.open(trajs.at(j).c_str(), std::ios::binary);

		if (trjin.is_open()){
			ftrjin=new Trajectory;
      ftrjin->setShow(show);

      if (pdb.length() > 0){
        ftrjin->setMolecule(mol);
      }

      if (ftrjin->findFormat(trjin) == true){
				ftrjin->readHeader(trjin);
				if (j == 0 && out.length() > 0 && trjout.is_open()){
					ftrjout->cloneHeader(ftrjin);
          ftrjout->writeHeader(trjout);
				}
				for (i=0; i< ftrjin->getNFrame(); i++){
					ftrjin->readFrame(trjin, i);
          if (pdb.length() >0){
            if (fit == true){

					  }
					  else if (recenter == true){

            }
            else{
              //Do nothing
            }
					}
          if (out.length() > 0 && trjout.is_open()){
            ftrjout->writeFrame(trjout, ftrjin, i);
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

  if (out.length() > 0 && trjout.is_open()){
    trjout.close();
		if (ftrjout != NULL){
			delete ftrjout;
		}
  }

  return 0;
}
