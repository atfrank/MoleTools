//Sean M. Law

#include "Residue.hpp"

Residue::Residue (){
  resname="   ";
  resid=0;
  chainid=" ";
  start=NULL;
  end=NULL;
  segid="    ";
  atmVec.clear();
}

void Residue::reset(){
  resname="   ";
  resid=0;
  chainid=" ";
  start=NULL;
  end=NULL;
  segid="    ";
  atmVec.clear();
}

void Residue::setResName(std::string resnamein){
  resname=resnamein;
}

void Residue::setResId(int residin){
  resid=residin;
}

void Residue::setChainId(std::string chainidin){
  chainid=chainidin;
}

void Residue::setStart(Atom* startin){
  start=startin;
}

void Residue::setEnd(Atom* endin){
  end=endin;
}

void Residue::setSegId(std::string segidin){
  segid=segidin;
}

void Residue::addAtom(Atom* atmEntry){
  if (atmEntry->getAtmNum()){
    atmVec.push_back(atmEntry);
  }
}

int Residue::getResId(){
  return resid;
}

std::string Residue::getResName(){
  return resname;
}

std::string Residue::getChainId(){
  return chainid;
}

std::string Residue::getSegId(){
  return segid;
}

Atom* Residue::getAtom (int element){
  return atmVec.at(element);
}

unsigned int Residue::getAtmVecSize (){
  return atmVec.size();
}
