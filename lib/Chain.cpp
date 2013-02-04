//Sean M. Law

#include "Chain.hpp"

Chain::Chain (){
  id.clear();
  resVec.clear();
  atmVec.clear();
}

void Chain::reset(){
  id.clear();
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

int Chain::getChainIdSize(){
  return id.size();
}

std::string Chain::getChainId(){
  return id;
}

void Chain::setChainId(const std::string& chainidin){
  id=chainidin;
}

void Chain::setChainId(){
  id.clear();
}

