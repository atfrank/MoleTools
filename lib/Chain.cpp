//Sean M. Law

#include "Chain.hpp"

Chain::Chain (){
  id.clear();
  atmVec.clear();
}

void Chain::reset(){
  id.clear();
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

