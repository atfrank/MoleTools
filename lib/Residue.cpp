//Sean M. Law

#include "Residue.hpp"

Residue::Residue (){
  resname=NULL;
  resid=NULL;
  chainid=NULL;
  start=NULL;
  end=NULL;
  segid=NULL;
  atmVec.clear();
  sel=true;
}

void Residue::reset(){
  resname=NULL;
  resid=NULL;
  chainid=NULL;
  start=NULL;
  end=NULL;
  segid=NULL;
  atmVec.clear();
  sel=true;
}

void Residue::addAtom(Atom* atmEntry){
  if (atmEntry->getAtmNum()){
    atmVec.push_back(atmEntry);
  }
}

int Residue::getResId(){
  return this->getAtom(0)->getResId();
}

std::string Residue::getResName(){
  return this->getAtom(0)->getResName();
}

std::string Residue::getChainId(){
  return this->getAtom(0)->getChainId();
}

Atom* Residue::getStart(){
  return this->getAtom(0);
}

Atom* Residue::getEnd(){
  return this->getAtom(atmVec.size()-1);
}

std::string Residue::getSegId(){
  return this->getAtom(0)->getSegId();
}

Atom* Residue::getAtom (int element){
  return atmVec.at(element);
}

unsigned int Residue::getAtmVecSize (){
  return atmVec.size();
}

void Residue::setSel (bool selin){
  sel=selin;
}

bool& Residue::getSel (){
  return sel;
}

void Residue::resetSel(){
  sel=true;
  for(unsigned int i=0; i< this->getAtmVecSize(); i++){
    this->getAtom(i)->setSel(true);
  }
}
