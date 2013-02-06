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

void Molecule::addResidue(Residue resEntry){
  if (resEntry.getResId()){
    resVec.push_back(resEntry);
  }
}

Chain* Molecule::getChain(int element){
  return &(chnVec.at(element));
}

unsigned int Molecule::getChnVecSize(){
  return chnVec.size();
}

Chain* Molecule::getLastChainRef(){
  if (chnVec.size() >0){
    return &(chnVec.at(chnVec.size()-1));
  }
  else{
    return NULL;
  }
}

unsigned int Molecule::getResVecSize(){
  return resVec.size();
}

Residue* Molecule::getLastResidueRef(){
  if (resVec.size() > 0){
    return &(resVec.at(resVec.size()-1));
  }
  else{
    return NULL;
  }
}

Residue* Molecule::getResidue(int element){
  return &(resVec.at(element));
}

Atom* Molecule::getLastAtomRef(){
  if (atmVec.size() >0){
    return &(atmVec.at(atmVec.size()-1));
  }
  else{
    return NULL;
  }
}
