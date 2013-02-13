//Sean M. Law

#include "Molecule.hpp"
#include "PDB.hpp"
#include "Select.hpp"

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

Molecule* Molecule::clone (bool selFlag){
  //Deep copy
  Molecule *mol=new Molecule;
  Chain *chnEntry=new Chain;
  Residue *resEntry=new Residue;
  Atom *atmEntry;
  Atom *lastAtom;

  atmEntry=NULL;
  lastAtom=NULL;

  //Create new Chains, Residues, Atoms
  for (unsigned int i=0; i< this->getAtmVecSize(); i++){
    if(selFlag == true && this->getAtom(i)->getSel() == false){
        continue;
    }
    atmEntry=new Atom; //This is necessary!
    atmEntry->clone(this->getAtom(i)); //Clone Atom

    /*****Same as PDB::readPDB*******/
    mol->addAtom(atmEntry);

    //Residue/Chain
    if (lastAtom == NULL){
      resEntry->addAtom(atmEntry);
      chnEntry->addAtom(atmEntry);
    }
    else if (atmEntry->getChainId() != lastAtom->getChainId()) {
      //Store last
      chnEntry->addResidue(resEntry);
      mol->addResidue(resEntry);
      mol->addChain(chnEntry);
      //Create new
      chnEntry=new Chain;
      chnEntry->addAtom(atmEntry);
      resEntry=new Residue;
      resEntry->addAtom(atmEntry);
    }
    else if (lastAtom->getResId() != atmEntry->getResId()){
      chnEntry->addResidue(resEntry);
      mol->addResidue(resEntry);
      resEntry=new Residue;
      resEntry->addAtom(atmEntry);
    }
    else{
      resEntry->addAtom(atmEntry);
    }

    //Update for next atom
    lastAtom=atmEntry;
  }
  chnEntry->addResidue(resEntry);
  mol->addResidue(resEntry);
  mol->addChain(chnEntry);
  /**************************/

  return mol;
}

Atom* Molecule::getAtom(int element){
  return atmVec.at(element);
}

unsigned int Molecule::getAtmVecSize(){
  return atmVec.size();
}

std::vector<Atom*> Molecule::getAtmVec(){
  return atmVec;
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

void Molecule::selAll(){
  unsigned int i;
//  for (i=0; i< this->getChnVecSize(); i++){
//    this->getChain(i)->setSel(true);
//  }
//  for (i=0; i< this->getResVecSize(); i++){
//    this->getResidue(i)->setSel(true);
//  }
  for (i=0; i< this->getAtmVecSize(); i++){
    this->getAtom(i)->setSel(true);
  }
}

void Molecule::deselAll(){
  unsigned int i;
//  for (i=0; i< this->getChnVecSize(); i++){
//    this->getChain(i)->setSel(false);
//  }
//  for (i=0; i< this->getResVecSize(); i++){
//    this->getResidue(i)->setSel(false);
//  }
  for (i=0; i< this->getAtmVecSize(); i++){
    this->getAtom(i)->setSel(false);
  }
}

void Molecule::select(std::string sel){
  Select::makeSel(this, sel);
}
