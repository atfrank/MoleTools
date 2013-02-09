//Sean M. Law

#include "Molecule.hpp"
#include "PDB.hpp"

#include <iostream>
#include <fstream>
#include <sstream> //string stream
#include <string>
#include <cstdlib>
using namespace std;

#define MAXINPFILES 4096
#define MAXLINESIZE 4096

void usage (){
  cerr << endl;
  cerr << "Usage:  manipPDB [options] <PDBfile>" << endl;
  cerr << "Options: [-model num]" << endl;
  cerr << endl << endl;
  exit(0);
}

int main (int argc, char **argv){

  int i;
  int model=0;
  string pdb;
  string currArg;

  pdb.clear();
  
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
      pdb=currArg;
    }
  }

  if (pdb.length() == 0){
    cerr << endl << "Error: Please provide an input file" << endl << endl;
    usage();
  }

  Molecule *mol=Molecule::readPDB(pdb, model);

  Molecule *cmol=mol->clone();
  cmol->resetSel();

  mol->writePDB();
//  cmol->writePDB();

  /*
  cout << mol->getAtom(0).getCoor().x() ;
  Atom t=mol->getAtom(0);
  cout << t.getCoor().x() << endl;
  cout << t.getX() << endl;
  cout << t.getAtmnum() << endl;
  */
  return 0;
}
