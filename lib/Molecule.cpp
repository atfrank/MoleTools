//Sean M. Law

#include "Molecule.hpp"
#include "PDB.hpp"

int Molecule::readPDB (std::string ifile, int model){
 
  PDB::readPDB(*this, ifile, model);

  return 0;
}

int Molecule::writePDB(){

  std::string out;

  out=PDB::writePDBFormat(*this);

  std::cout << out;

  return 0;
}

void Molecule::addAtom(Atom atmEntry) {
  if (atmEntry.getAtmNum()){
    atmVec.push_back(atmEntry);
  }
}

Atom Molecule::getAtom(int atmnum){
  return atmVec.at(atmnum);
}

unsigned int Molecule::getAtmVecSize(){
  return atmVec.size();
}

void Molecule::addChain(Chain chnEntry){
  chnVec.push_back(chnEntry);
}
