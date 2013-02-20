//Sean M. Law

#include "SABA.hpp"

Molecule* SABA::getPseudoCenter(Molecule *mol){
	Molecule *ssmol=new Molecule;
	Chain *chnEntry=new Chain;
	Residue *resEntry=new Residue;
	Atom *atmEntry;
	Atom *lastAtom;
  Vector xyz;

	atmEntry=NULL;
	lastAtom=NULL;

	for (unsigned int i=0; i< mol->getChnVecSize(); i++){
		for (unsigned int j=0; j< mol->getChain(i)->getAtmVecSize(); j++){

			atmEntry=new Atom;
			atmEntry->clone(mol->getChain(i)->getAtom(j));

			/*****Same as PDB::readPDB*******/
			if (j < mol->getAtmVecSize() - 1){
    		ssmol->addAtom(atmEntry);
        if (lastAtom != NULL){
          //Pseudo center
          //xyz=(lastAtom->getCoor() + atmEntry->getCoor())/2.0;
          //lastAtom->setCoor(xyz);
        }
			}
			else{
				//Clean up last by getting pseudo center
        //xyz=(lastAtom->getCoor() + atmEntry->getCoor())/2.0;
        //lastAtom->setCoor(xyz);
        continue;
			}

    	//Residue/Chain
    	if (lastAtom != NULL && atmEntry->getChainId() != lastAtom->getChainId()) {
      	//Store last
      	chnEntry->addResidue(resEntry);
      	ssmol->addResidue(resEntry);
      	ssmol->addChain(chnEntry);
      	//Create new
      	chnEntry=new Chain;
      	resEntry=new Residue;
    	}
    	else if (lastAtom != NULL && lastAtom->getResId() != atmEntry->getResId()){
      	chnEntry->addResidue(resEntry);
      	ssmol->addResidue(resEntry);
      	resEntry=new Residue;
    	}
    	else{
				//Do nothing
    	}

      if (j < mol->getAtmVecSize() - 1){
			  chnEntry->addAtom(atmEntry);
    	  resEntry->addAtom(atmEntry);
      }
			
    	//Update for next atom
    	lastAtom=atmEntry;		
		}
	}
	chnEntry->addResidue(resEntry);
  ssmol->addResidue(resEntry);
  ssmol->addChain(chnEntry);
  /**************************/

	return ssmol;
}
