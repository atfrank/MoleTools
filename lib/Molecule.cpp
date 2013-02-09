//Sean M. Law

#include "Molecule.hpp"
#include "PDB.hpp"

Molecule* Molecule::readPDB (std::string ifile, int model){

  return PDB::readPDB(ifile, model);

}

int Molecule::writePDB(){

  std::string out;

  out=PDB::writePDBFormat(this);

  std::cout << out;

  return 0;
}

void Molecule::addAtom(Atom* atmEntry) {
  if (atmEntry->getAtmNum()){
    atmVec.push_back(atmEntry);
  }
}

Molecule* Molecule::clone (){
  //Deep copy
  Molecule *mol=new Molecule;
  Atom *lastAtom;

  lastAtom=NULL;

  //Create new Chains, Residues, Atoms
  for (unsigned int i=0; i< this->getAtmVecSize(); i++){
    Atom *atmEntry=new Atom;
    atmEntry->clone(this->getAtom(i)); //Clone Atom
    mol->addAtom(atmEntry);

  }

  //Adjust pointers
  return mol;
}

Atom* Molecule::getAtom(int element){
  return atmVec.at(element);
}

unsigned int Molecule::getAtmVecSize(){
  return atmVec.size();
}

void Molecule::addChain(Chain* chnEntry){
  if(chnEntry->getChainId().size()){ //Check if empty string
    chnVec.push_back(chnEntry);
  }
}

void Molecule::addResidue(Residue* resEntry){
  if (resEntry->getResId()){
    resVec.push_back(resEntry);
  }
}

Chain* Molecule::getChain(int element){
  return chnVec.at(element);
}

unsigned int Molecule::getChnVecSize(){
  return chnVec.size();
}

unsigned int Molecule::getResVecSize(){
  return resVec.size();
}

Residue* Molecule::getResidue(int element){
  return resVec.at(element);
}
