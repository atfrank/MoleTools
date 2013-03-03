//Sean M. Law

#include "Trajectory.hpp"
#include "Molecule.hpp"

#include <iostream>
#include <fstream>
#include <cstdlib>
using namespace std;

#define MAXINPFILES 4096
#define MAXLINESIZE 4096

void usage(){
  cerr << endl << endl;
  cerr << "Usage:   modTraj [-options] <TRAJfile(s)>" << endl;
  cerr << "Options: [-pdb PDBfile] [-sel selection]" << endl;
  cerr << "         [-fit | -recenter] [-fitsel selection] [-recsel selection]" << endl;
  cerr << "         [-out TRAJname] [-outsel selection]" << endl;
  cerr << "         [-skip frames]" << endl;
	cerr << "         [-verbose] [-show]" << endl;
  exit(0);
}

int main (int argc, char **argv){


  int i;
	unsigned int j;
  vector<string> trajs;
  string out;
  Molecule *mol;
  Molecule *outmol;
  Molecule *fitmol;
  Molecule *recmol;
  string pdb;
  string currArg;
  string sel;
  bool fit=false;
  string fitsel=":.CA+P";
  bool recenter=false;
  string recsel=":.CA+P";
  string outsel=":.";
	ifstream trjin;
	ofstream trjout;
	Trajectory *ftrjin;
	Trajectory *ftrjout;

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
    else if (currArg == "-sel" || currArg == "-nsel"){
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
    else{
      trajs.push_back(currArg);
    }
  }

  if (trajs.size() == 0){
    cerr << endl << "Error: Please provide an input trajectory file" << endl << endl;
    usage();
  }
  else if (outsel.length() > 0 && out.length() == 0){
    cerr << endl << "Error: Please specify an output trajectory via \"-out\" ";
    cerr << endl << "when using option \"-outsel\"" << endl;
  }

  if (fit == true && recenter == true){
    cerr << endl << "Error: Options \"-fit\" and \"-recenter\" cannot be used simultaneously";
    cerr << endl;
    usage();
  }

  if (sel.length() > 0 && pdb.length() == 0){
    cerr << endl << "Error: Please provide a PDB file via \"-pdb\"";
    cerr << "in order to make a selection" << endl;
    usage();
  }
  else if (pdb.length() > 0){
    mol=Molecule::readPDB(pdb);

    if (outsel.length() > 0){
      outmol=mol->clone();
      //outmol->select(outsel);
    }
    if (fitsel.length() > 0 && fit == true){
      fitmol=mol->clone();
      fitmol->select(fitsel);
    }
    if (recsel.length() > 0 && recenter == true){
      recmol=mol->clone();
      recmol->select(recsel);
    }
    if (sel.length() > 0){
      mol->select(sel); //Currently serves no function
    }
  } 

  if (out.length() > 0){
    trjout.open(out.c_str(), ios::binary);
		ftrjout=new Trajectory;
    if (outmol != NULL){
      ftrjout->setMolecule(outmol);
    }
  }

	for (j=0; j< trajs.size(); j++){
		trjin.open(trajs.at(j).c_str(), ios::binary);

		if (trjin.is_open()){
			ftrjin=new Trajectory;

      if (pdb.length() > 0){
        //ftrjin->setMolecule(mol);
      }

      if (ftrjin->findFormat(trjin) == true){
				ftrjin->readHeader(trjin);
				//ftrjin->showHeader();
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
				cerr << "Warning: Skipping unknown trajectory format \"";
				cerr << trajs.at(j) << "\"" << endl;
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
