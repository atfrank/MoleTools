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

Atom* Residue::getAtom (int element){
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

	if(aa == "ALA"){return "A";}
	else if(aa ==	"CYS"){return "C";}
	else if (aa == "ASP"){return "D";}
	else if (aa == "GLU"){return "E";}
	else if (aa == "PHE"){return "F";}
	else if (aa == "GLY"){return "G";}
	else if (aa == "HIS" || aa == "HSD" || aa == "HSE" || aa == "HSP"){
		return "H";
	}
	else if (aa == "ILE"){return "I";}
	else if (aa == "LYS"){return "K";}
	else if (aa == "LEU"){return "L";}
	else if (aa == "MET"){return "M";}
	else if (aa == "ASN"){return "N";}
	else if (aa == "PRO"){return "P";}
	else if (aa == "GLN"){return "Q";}
	else if (aa == "ARG"){return "R";}
	else if (aa == "SER"){return "S";}
	else if (aa == "THR"){return "T";}
	else if (aa == "VAL"){return "V";}
	else if (aa == "TRP"){return "W";}
	else if (aa == "TYR"){return "Y";}
	else{
		return "";
	}
}

std::string Residue::aa123(const std::string &aa){
	if(aa == "A"){return "A";}
  else if(aa == "C"){return "CYS";} 
  else if (aa == "D"){return "ASP";}
  else if (aa == "E"){return "GLY";}
  else if (aa == "F"){return "PHE";}
  else if (aa == "G"){return "GLY";}
  else if (aa == "H"){return "HIS";}
  else if (aa == "I"){return "ILE";}
  else if (aa == "K"){return "LYS";}
  else if (aa == "L"){return "LEU";}
  else if (aa == "M"){return "MET";}
  else if (aa == "N"){return "ASN";}
  else if (aa == "P"){return "PRO";}
  else if (aa == "Q"){return "GLN";}
  else if (aa == "R"){return "ARG";}
  else if (aa == "S"){return "SER";}
  else if (aa == "T"){return "THR";}
  else if (aa == "V"){return "VAL";}
  else if (aa == "W"){return "TRP";}
  else if (aa == "Y"){return "TYR";}
  else{
    return "";
  }	
}
