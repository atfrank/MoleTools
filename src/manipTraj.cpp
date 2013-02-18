//Sean M. Law

#include "Molecule.hpp"

#include <iostream>
#include <fstream>
#include <cstdlib>
using namespace std;

#define MAXINPFILES 4096
#define MAXLINESIZE 4096

void usage(){
  
  exit(0);
}

int main (int argc, char **argv){


  int i;
  vector<string> dcds;
  string pdb;
  string currArg;
  string sel;

  pdb.clear();

  for (i=1; i<argc; i++){
    currArg=argv[i];
    if (currArg == "-h" || currArg == "-help"){
      usage();
    }
    else if (currArg == "-nsel"){
      currArg=argv[++i];
      sel=currArg;
    }
    else if (currArg == "-pdb"){
      currArg=argv[++i];
      pdb=currArg;
    }
    else{
      dcds.push_back(currArg);
    }
  }

  if (dcds.size() == 0){
    cerr << endl << "Error: Please provide an input trajectory file" << endl << endl;
    usage();
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

  return 0;
}
