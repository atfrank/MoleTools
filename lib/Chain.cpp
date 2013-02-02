//Sean M. Law

#include "Chain.hpp"

Chain::Chain (){
  id.clear();
  atmVec.clear();
}

void Chain::reset(){
  id.clear();
}

void Chain::setChainId(const std::string& chainidin){
  id=chainidin;
}

void Chain::setChainId(){
  id.clear();
}

