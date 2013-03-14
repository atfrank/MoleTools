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
  std::cerr << "Usage:   modTraj [-options] <TRAJfile(s)>" << std::endl;
  std::cerr << "Options: [-pdb PDBfile] [-sel selection]" << std::endl;
  std::cerr << "         [-fit refPDB] [-fitsel selection] [-recsel selection]" << std::endl;
  std::cerr << "         [-out TRAJname] [-outsel selection]" << std::endl;
  std::cerr << "         [-skip frames] [-start frame]" << std::endl;
	std::cerr << "         [-center]" << std::endl;
	std::cerr << "         [-translate dx dy dz]" << std::endl;
	std::cerr << "         [-show]" << std::endl;
  exit(0);
}

int main (int argc, char **argv){

  int i;
	unsigned int j;
  std::vector<std::string> trajs;
  std::string fout;
	bool out=false;
  Molecule *mol;
  Molecule *outmol;
  Molecule *refmol;
  Molecule *recmol;
  std::string pdb;
	std::string refpdb;
  std::string sel=""; //To match NATOM in trajectory
  std::string currArg;
  bool fit=false;
  std::string fitsel=":.";
  bool recenter=false;
  std::string recsel=":.";
	bool translate=false;
	double dx, dy, dz;
	Vector dxyz;
	bool center=false;
  std::string outsel="";
  std::ifstream trjin;
  std::ofstream trjout;
	Trajectory *ftrjin;
	Trajectory *ftrjout;
  bool show=false;
  int skip=0;
  int start=0;

  pdb.clear();
	refpdb.clear();
  mol=NULL;
  outmol=NULL;
  refmol=NULL;
  recmol=NULL;
	ftrjin=NULL;
	ftrjout=NULL;
	dx=0;
	dy=0;
	dz=0;

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
			currArg=argv[++i];
			refpdb=currArg;
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
		else if (currArg == "-center"){
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
    else if (currArg == "-out"){
      currArg=argv[++i];
      fout=currArg;
			out=true;
    }
    else if (currArg == "-outsel"){
      currArg=argv[++i];
      outsel=currArg;
			out=true;
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
  else if (out == true && fout.length() == 0){
    std::cerr << std::endl << "Error: Please specify an output trajectory via \"-out\" ";
    std::cerr << std::endl << "when using option \"-outsel\"" << std::endl;
  }

  if (fit == true && recenter == true){
    std::cerr << std::endl << "Error: Options \"-fit\" and \"-recenter\" cannot be used simultaneously";
    std::cerr << std::endl;
    usage();
  }

  if (((out == true && outsel.length() > 0)) && pdb.length() == 0){
    std::cerr << std::endl << "Error: Please provide a PDB file via \"-pdb\"" << std::endl;
    usage();
  }
	else if (pdb.length() == 0 && (recenter == true || fit == true)){
		std::cerr << std::endl << "Error: Please provide a PDB file via \"-pdb\"" << std::endl;
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
			if (refpdb.length() == 0){
				std::cerr << std::endl << "Error: Please provide a reference PDB file via \"-fit\"" << std::endl;
    		usage();
			}
			refmol=Molecule::readPDB(refpdb);
			refmol->select(fitsel);
      refmol=refmol->clone(true, false); //Delete original after cloning
    }
    if (recenter == true){
      mol->select(recsel);
      recmol=mol->copy();
    }
  }
  else{
    //Do Nothing
  }

  if (out == true){
    trjout.open(fout.c_str(), std::ios::binary);
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
				if (fit == true){
					ftrjin->getMolecule()->select(fitsel);
				}
      }

      if (ftrjin->findFormat(trjin) == true){
				ftrjin->readHeader(trjin);
				if (j == 0 && out == true && trjout.is_open()){
					ftrjout->cloneHeader(ftrjin);
          ftrjout->setNFrame(0);
          if (start != 0){
            ftrjout->setNPriv(ftrjout->getNPriv()+start*ftrjout->getNSavc());
          }
          if (skip != 0){
            ftrjout->setNSavc(ftrjout->getNSavc()*(skip+1));
          }
          ftrjout->writeHeader(trjout);
				}
        //Loop through desired frames
				for (i=start; i< ftrjin->getNFrame(); i=i+1+skip){
					ftrjin->readFrame(trjin, i);
          if (pdb.length() >0){
						//Transform molecule only if PDB is provided
            if (fit == true){
							ftrjin->getMolecule()->lsqfit(refmol);
					  }
					  else if (recenter == true){
							ftrjin->getMolecule()->recenter(recmol);
            }
						else if (center == true){
							ftrjin->getMolecule()->center();
						}
						else if (translate == true){
							ftrjin->getMolecule()->translate(dxyz);
						}
            else{
              //Do nothing
            }
					}
          if (out == true && trjout.is_open()){
            ftrjout->setNFrame(ftrjout->getNFrame()+1);
            ftrjout->writeFrame(trjout, ftrjin);
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

  if (out == true && trjout.is_open()){
    //Re-write header in case of any changes
    ftrjout->writeHeader(trjout);
    trjout.close();
		if (ftrjout != NULL){
			delete ftrjout;
		}
  }

  return 0;
}
