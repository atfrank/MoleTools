//Sean M. Law

#include "Molecule.hpp"
#include "PDB.hpp"
#include "Select.hpp"

Molecule::~Molecule (){
	Chain *c;
	Residue *r;
	Atom *a;

	for (unsigned int i=0; i< this->getChnVecSize(); i++){
		c=this->getChain(i);
		for (unsigned int j=0; j< c->getResVecSize(); j++){
			r=c->getResidue(j);
			for (unsigned int k=0; k< r->getAtmVecSize(); k++){
				a=r->getAtom(k);
				delete a;
			}
			delete r;
		}
		delete c;
	}
}

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
  Chain *c, *chnEntry;
  Residue *r, *resEntry;
  Atom *a, *atmEntry;

	c=NULL;
	r=NULL;
	a=NULL;
	chnEntry=NULL;
	resEntry=NULL;
  atmEntry=NULL;
	
  for (unsigned int i=0; i< this->getChnVecSize(); i++){
    c=this->getChain(i);
    chnEntry=new Chain; //Deleted later if no residues/atoms

    for (unsigned int j=0; j< c->getResVecSize(); j++){
      r=c->getResidue(j);
      resEntry=new Residue; //Deleted later if no atoms

      for (unsigned int k=0; k< r->getAtmVecSize(); k++){
        a=r->getAtom(k);
        if(selFlag == true && a->getSel() == false){
          continue;
        }
        //Add each selected atom
				atmEntry=new Atom; //This is necessary!
				atmEntry->clone(a); //Clone Atom
        mol->addAtom(atmEntry);
        resEntry->addAtom(atmEntry);
        chnEntry->addAtom(atmEntry);
      }
      if (resEntry->getAtmVecSize() > 0){
        //Add each residue with selected atoms
        mol->addResidue(resEntry);
        chnEntry->addResidue(resEntry);
      }
      else{
        delete resEntry;
      }
    }
    if (chnEntry->getResVecSize() > 0){
      //Add each chain with selected atoms
      mol->addChain(chnEntry);
    }
    else{
      delete chnEntry;
    }
  }

  return mol;
}

Molecule* Molecule::copy (bool selFlag){
  //Shallow copy of selected atom pointers.
	//Chains and Residues are required to be created new.
	//Useful when you need to make multiple selections
  Molecule *mol=new Molecule;
  Chain *c, *chnEntry;
  Residue *r, *resEntry;
  Atom *a;

	c=NULL;
	r=NULL;
  a=NULL;
	chnEntry=NULL;
	resEntry=NULL;

	for (unsigned int i=0; i< this->getChnVecSize(); i++){
		c=this->getChain(i);
		chnEntry=new Chain;

		for (unsigned int j=0; j< c->getResVecSize(); j++){
			r=c->getResidue(j);
			resEntry=new Residue;

			for (unsigned int k=0; k< r->getAtmVecSize(); k++){
				a=r->getAtom(k);
				if(selFlag == true && a->getSel() == false){
     	 		continue;
    		}
				//Add each selected atom
				mol->addAtom(a); //Copy atom pointer, not clone
				resEntry->addAtom(a);
				chnEntry->addAtom(a);
			}
			if (resEntry->getAtmVecSize() > 0){
				//Add each residue with selected atoms
				mol->addResidue(resEntry);
				chnEntry->addResidue(resEntry);
			}
			else{
				delete resEntry;
			}
		}
		if (chnEntry->getResVecSize() > 0){
			//Add each chain with selected atoms
			mol->addChain(chnEntry);
		}
		else{
			delete chnEntry;
		}			
	}	

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

unsigned int Molecule::getNatom(){
  return atmVec.size();
}
