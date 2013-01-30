#include <iostream>
#include <fstream>
#include <sstream> //string stream
#include <string>
#include <cstdlib>
using namespace std;

#include "Molecule.hpp"

#define MAXINPFILES 4096
#define MAXLINESIZE 4096

void usage (){
  cerr << endl;
  cerr << "Usage:  manipPDB [options] <PDBfile>" << endl;
  cerr << "Options: [-h || -help]" << endl;
  cerr << endl << endl;
  exit(0);
}

int main (int argc, char **argv){

  int i;
  int model;
  string *pdb = new string;
  string currArg;
  
  for (i=1; i<argc; i++){
    currArg=argv[i];
    if (currArg == "-h" || currArg == "-help"){
      usage();
    }
    else if (currArg == "-model"){
      currArg=argv[++i];
      stringstream(currArg) >> model; //atoi
    }
    else{
      *pdb=currArg;
    }
  }

//  Molecule *mol = new Molecule;
//  (*mol).readPDB(pdb);

  return 0;
}
