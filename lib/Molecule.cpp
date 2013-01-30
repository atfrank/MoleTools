#include "Molecule.hpp"

Molecule::Molecule(){
  
}

Molecule::Molecule(int atmnum, char *atmname, char *resname, int resnum, Vector vec, char *seg){

}

void Molecule::readPDB (string *ifile){
  
  ifstream PDBfile;

  if (*ifile == "-"){

  }
  else{
    PDBfile.open((*ifile).c_str());
  }
}
