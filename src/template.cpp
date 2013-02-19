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
    else{
      pdb=currArg;
    }
  }

  if (pdb.length() == 0){
    cerr << endl << "Error: Please provide an input file" << endl << endl;
    usage();
  }

  Molecule *mol=Molecule::readPDB(pdb);

  if (sel.length() >0){
    mol->select(sel);
  }

  Molecule *cmol=mol->clone();
  

  return 0;
}
