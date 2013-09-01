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
//  sel=true;
}

void Residue::reset(){
  resname=NULL;
  resid=NULL;
  chainid=NULL;
  start=NULL;
  end=NULL;
  segid=NULL;
  atmVec.clear();
//  sel=true;
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

Atom* Residue::getAtom (const int &element){
  return atmVec.at(element);
}

unsigned int Residue::getAtmVecSize (){
  return atmVec.size();
}

//void Residue::setSel (bool selin){
//  sel=selin;
//}

//bool& Residue::getSel (){
//  return sel;
//}

void Residue::selAll(){
//  sel=true;
  for(unsigned int i=0; i< this->getAtmVecSize(); i++){
    this->getAtom(i)->setSel(true);
  }
}

void Residue::deselAll(){
//  sel=false;
  for(unsigned int i=0; i< this->getAtmVecSize(); i++){
    this->getAtom(i)->setSel(false);
  }
}

std::string Residue::aa321(const std::string &aa){

	if(aa.compare(0,3,"ALA") == 0){return "A";}
	else if(aa.compare(0,3,"CYS") == 0){return "C";}
	else if (aa.compare(0,3,"ASP") == 0){return "D";}
	else if (aa.compare(0,3,"GLU") == 0){return "E";}
	else if (aa.compare(0,3,"PHE") == 0){return "F";}
	else if (aa.compare(0,3,"GLY") == 0){return "G";}
	else if (aa.compare(0,3,"HIS") == 0 || aa.compare(0,3,"HSD") == 0 || aa.compare(0,3,"HSE") == 0 || aa.compare(0,3,"HSP") == 0){
		return "H";
	}
	else if (aa.compare(0,3,"ILE") == 0){return "I";}
	else if (aa.compare(0,3,"LYS") == 0){return "K";}
	else if (aa.compare(0,3,"LEU") == 0){return "L";}
	else if (aa.compare(0,3,"MET") == 0){return "M";}
	else if (aa.compare(0,3,"ASN") == 0){return "N";}
	else if (aa.compare(0,3,"PRO") == 0){return "P";}
	else if (aa.compare(0,3,"GLN") == 0){return "Q";}
	else if (aa.compare(0,3,"ARG") == 0){return "R";}
	else if (aa.compare(0,3,"SER") == 0){return "S";}
	else if (aa.compare(0,3,"THR") == 0){return "T";}
	else if (aa.compare(0,3,"VAL") == 0){return "V";}
	else if (aa.compare(0,3,"TRP") == 0){return "W";}
	else if (aa.compare(0,3,"TYR") == 0){return "Y";}
	else{
		return "";
	}
}

std::string Residue::aa123(const std::string &aa){
	if(aa.compare(0,1,"A") == 0){return "A";}
  else if(aa.compare(0,1,"C") == 0){return "CYS";} 
  else if (aa.compare(0,1,"D") == 0){return "ASP";}
  else if (aa.compare(0,1,"E") == 0){return "GLY";}
  else if (aa.compare(0,1,"F") == 0){return "PHE";}
  else if (aa.compare(0,1,"G") == 0){return "GLY";}
  else if (aa.compare(0,1,"H") == 0){return "HIS";}
  else if (aa.compare(0,1,"I") == 0){return "ILE";}
  else if (aa.compare(0,1,"K") == 0){return "LYS";}
  else if (aa.compare(0,1,"L") == 0){return "LEU";}
  else if (aa.compare(0,1,"M") == 0){return "MET";}
  else if (aa.compare(0,1,"N") == 0){return "ASN";}
  else if (aa.compare(0,1,"P") == 0){return "PRO";}
  else if (aa.compare(0,1,"Q") == 0){return "GLN";}
  else if (aa.compare(0,1,"R") == 0){return "ARG";}
  else if (aa.compare(0,1,"S") == 0){return "SER";}
  else if (aa.compare(0,1,"T") == 0){return "THR";}
  else if (aa.compare(0,1,"V") == 0){return "VAL";}
  else if (aa.compare(0,1,"W") == 0){return "TRP";}
  else if (aa.compare(0,1,"Y") == 0){return "TYR";}
  else{
    return "";
  }	
}
