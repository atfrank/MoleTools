//Sean M. Law

#include "Molecule.hpp"
#include "SABA.hpp"

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

  pdb.clear();

  for (i=1; i<argc; i++){
    currArg=argv[i];
    if (currArg == "-h" || currArg == "-help"){
      usage();
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

  mol->select(":.CA");
  mol=mol->clone();
  mol->writePDB(); 

  return 0;
}
