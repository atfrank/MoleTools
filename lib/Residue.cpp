//Sean M. Law

#include "Residue.hpp"

Residue::Residue (){
  resname="   ";
  resnum=0;
  chainid=" ";
  start=NULL;
  end=NULL;
  segid="    ";
  atmVec.clear();
}

void Residue::setResName(std::string resnamein){
  resname=resnamein;
}

void Residue::setResNum(int resnumin){
  resnum=resnumin;
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

int Residue::getResNum(){
  return resnum;
}
