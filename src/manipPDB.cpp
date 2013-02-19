//Sean M. Law

#include "Molecule.hpp"

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
  cerr << "         [-sel selection]" << endl;
  cerr << endl << endl;
  exit(0);
}

int main (int argc, char **argv){

  int i;
  int model=0;
  string pdb;
  string currArg;
  string sel;

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
    else if (currArg == "-sel" || currArg == "-nsel"){
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

  Molecule *mol=Molecule::readPDB(pdb, model);
  Molecule *cmol=mol->clone();
	
  cmol->deselAll();
  cmol->selAll();

  if (sel.length() >0){
    cmol->select(sel);
  }

//  mol->writePDB();
  cmol->writePDB();

  /*
  cout << mol->getAtom(0).getCoor().x() ;
  Atom t=mol->getAtom(0);
  cout << t.getCoor().x() << endl;
  cout << t.getX() << endl;
  cout << t.getAtmnum() << endl;
  */
  return 0;
}
