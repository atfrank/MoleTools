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

Atom Molecule::getAtom(int element){
  return atmVec.at(element);
}

unsigned int Molecule::getAtmVecSize(){
  return atmVec.size();
}

void Molecule::addChain(Chain chnEntry){
  if(chnEntry.getChainIdSize()){
    chnVec.push_back(chnEntry);
  }
}

Chain Molecule::getChain(int element){
  return chnVec.at(element);
}

unsigned int Molecule::getChnVecSize(){
  return chnVec.size();
}

Atom& Molecule::getLastAtomRef(){
  return atmVec.at(atmVec.size()-1);
}
