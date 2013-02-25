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
  cerr << "         [-fit] [-fitsel selection]" << endl;
  cerr << "         [-recenter] [-recsel selection]" << endl;
  cerr << "         [-out TRAJname] [-outsel selection]" << endl;
  exit(0);
}

int main (int argc, char **argv){


  int i;
	unsigned int j;
  vector<string> trajs;
  string out;
  string pdb;
  string currArg;
  string sel;
  bool fit=false;
  string fitsel=":.CA+P";
  bool recenter=false;
  string recsel=":.CA+P";
  string outsel="";
	ifstream trjin;
	ofstream trjout;

  pdb.clear();

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

  if (sel.length() > 0 && pdb.length()){
    cerr << endl << "Error: Please provide a PDB file via \"-pdb\"";
    cerr << "in order to make a selection" << endl;
    usage();
  }
  else{
    Molecule *mol=Molecule::readPDB(pdb);
    mol->select(sel); //Don't need to clone, just write selected
  } 

  if (out.length() > 0){
    trjout.open(out.c_str(), ios::binary);
  }

	for (j=0; j< trajs.size(); j++){
		trjin.open(trajs.at(j).c_str(), ios::binary);

		if (trjin.is_open()){
			
		}
	
		trjin.close();
	}

  if (out.length() > 0){
    trjout.close();
  }

  return 0;
}
