//Sean M. Law

#include "Chain.hpp"

Chain::Chain (){
  id=NULL;
  resVec.clear();
  atmVec.clear();
}

void Chain::reset(){
  id=NULL;
  resVec.clear();
  atmVec.clear();
}

void Chain::addResidue(Residue* resEntry){
  if (resEntry->getResId()){
    resVec.push_back(resEntry);
  }
}

void Chain::addAtom(Atom* atmEntry){
  if (atmEntry->getAtmNum()){
    atmVec.push_back(atmEntry);
  }
}

Atom* Chain::getAtom (int element){
  return atmVec.at(element);
}

Residue* Chain::getResidue (int element){
  return resVec.at(element);
}

std::string Chain::getChainId(){
  return this->getAtom(0)->getChainId();
}

unsigned int Chain::getAtmVecSize(){
  return atmVec.size();
}

unsigned int Chain::getResVecSize(){
  return resVec.size();
}
