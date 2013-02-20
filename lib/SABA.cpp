//Sean M. Law

#include "SABA.hpp"

Molecule* SABA::getPseudoCenter(Molecule *mol){
	Molecule *ssmol=new Molecule;
	Atom *atmEntry;
	Atom *lastAtom;
	Vector xyz;
	unsigned int i, j, k;

	//Leave original molecule untouched
	ssmol=mol->clone();
	ssmol->select(":.CA");
	ssmol=ssmol->clone();

	//Convert C-alpha to pseudo centers, calculate C-alpha distances
	for (i=0; i< ssmol->getChnVecSize(); i++){
		//Important to nullify for each chain
		atmEntry=NULL;
		lastAtom=NULL;
	  for (j=0; j< ssmol->getChain(i)->getAtmVecSize(); j++){
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
		//No i+1 neighbor, modify coordinates for last atom in chain
		xyz=Vector(9999.0, 9999.0, 9999.0);
		lastAtom->setCoor(xyz);
	}


	Atom *atm1, *atm2, *atm3, *atm4;
	Chain *chn1;
	int atmVecSize;
	std::vector<double>::iterator it;

	//Calculate distances and dihedral angles
	for (i=0; i< ssmol->getChnVecSize(); i++){
		chn1=ssmol->getChain(i);
		atmVecSize=chn1->getAtmVecSize();
		for (j=0; j< atmVecSize; j++){
			atm1=chn1->getAtom(j);
			atm2=NULL;
			atm3=NULL;
			atm4=NULL;
			if (j+1 < atmVecSize && atm1->getResId()+1 == chn1->getAtom(j+1)->getResId()){
				atm2=chn1->getAtom(j+1);
			}
			if (j+2 < atmVecSize && atm1->getResId()+2 == chn1->getAtom(j+2)->getResId()){
				atm3=chn1->getAtom(j+2);
			}
			if (j+3 < atmVecSize && atm1->getResId()+3 == chn1->getAtom(j+3)->getResId()){
				atm4=chn1->getAtom(j+3);
			}
			//Perform calculations
			if (atm1 != NULL && atm2 != NULL && atm3 != NULL && atm4 != NULL){
				
			}
		}
	}


	return ssmol;
}
