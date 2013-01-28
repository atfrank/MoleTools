#include <iostream>
#include <fstream>
#include <sstream> //string stream
#include <string>
using namespace std;

#include <Molecule.hpp>
#include <PDB.hpp>

#define MAXINPFILES 4096
#define MAXLINESIZE 4096

void usage (){
  cerr << endl;
  cerr << "Usage:  manipPDB [options] <PDBfile>" << endl;
  cerr << "Options: [-h || -help]" << endl;
  cerr << endl << endl;
  return;
}

int main (int argc, char **argv){

  ifstream fpdb;
  int i;
  int model;
  string pdb;
  string currArg;
  stringstream ss;
  
  for (i=1; i<argc; i++){
    currArg=argv[i];
    if (currArg == "-h" || currArg == "-help"){
      usage();
    }
    else if (currArg == "-model"){
      currArg=argv[++i];
      ss << currArg;
      ss >> model; //atoi
    }
    else{
      pdb=currArg;
    }
  }
  
//  fpdb.open(pdb.c_str());

  if (fpdb){
//    cerr << "File: " << pdb;
  }

  return 0;
}
