//Sean M. Law

#include "SABA.hpp"

Molecule* SABA::getPseudoCenter(Molecule *mol){
	Molecule *ssmol=new Molecule;
	Atom *atmEntry;
	Atom *lastAtom;
	Vector xyz;

	//Leave original molecule untouched
	ssmol=mol->clone();
	ssmol->select(":.CA");
	ssmol=ssmol->clone();

	for (unsigned int i=0; i< ssmol->getChnVecSize(); i++){
		//Important to nullify for each chain
		atmEntry=NULL;
		lastAtom=NULL;
	  for (unsigned int j=0; j< ssmol->getChain(i)->getAtmVecSize(); j++){
			atmEntry=ssmol->getChain(i)->getAtom(j);

			if (lastAtom != NULL){
				if(lastAtom->getResId()+1 == atmEntry->getResId()){
					//Computer pseudo center coordinates
					xyz=(lastAtom->getCoor() + atmEntry->getCoor())/2.0;
					lastAtom->setCoor(xyz);
				}
				else{
					//No i+1 neighbor, modify coordinates
					xyz=Vector(9999.0, 9999.0, 9999.0);
				}
				lastAtom->setCoor(xyz);
			}

			lastAtom=atmEntry;
		}
	}

	return ssmol;
}
