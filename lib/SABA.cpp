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
	unsigned int size;
	double dist1, dist2, angl, dihe;

  std::vector<double> helixDist; //Alpha/310 - i', i'+3
  std::vector<double> helixDihe; //Alpha/310 - i', i'+1, i'+2, i'+3
  std::vector<double> threeDist1; //310 - i'+1, i'+2
  std::vector<double> threeDist2; //310 - i'+1, i'+3
  std::vector<double> betaDist; //Parallel/Anti - i', j'
  std::vector<double> paraDist; //Parallel - i'-1, j'-1
  std::vector<double> antiDist1; //Anti - i'+1, j'-1
  std::vector<double> antiDist2; //Anti - i+1, j (C-alpha)

	//Calculate distances and dihedral angles
	for (i=0; i< ssmol->getChnVecSize(); i++){
		chn1=ssmol->getChain(i);
		size=chn1->getAtmVecSize();

		for (j=0; j< size; j++){
			atm1=chn1->getAtom(j);
			atm2=NULL;
			atm3=NULL;
			atm4=NULL;
			if (j+1 < size && atm1->getResId()+1 == chn1->getAtom(j+1)->getResId()){
				atm2=chn1->getAtom(j+1);
			}
			if (j+2 < size && atm1->getResId()+2 == chn1->getAtom(j+2)->getResId()){
				atm3=chn1->getAtom(j+2);
			}
			if (j+3 < size && atm1->getResId()+3 == chn1->getAtom(j+3)->getResId()){
				atm4=chn1->getAtom(j+3);
			}
			//Perform calculations
			dist1=999.0;
			dist2=999.0;
			angl=999.0;
			dihe=999.0;
			if (atm1 != NULL && atm2 != NULL && atm3 != NULL && atm4 != NULL){
				dist1=Vector::distance(atm1->getCoor(), atm4->getCoor());
				dist2=Vector::distance(atm1->getCoor(), atm3->getCoor());
				
				if (atm1->getResName() == "PRO"){
					//Do nothing
				}	
				else if (dist1 > 4.21 && dist1 < 5.23){
					dihe=Vector::dihedral(atm1->getCoor(), atm2->getCoor(),atm3->getCoor(),atm4->getCoor());
					if (dihe > 43.5 && dihe < 78.3){
						//Alpha Helix
						atm1->setSS("H");
					}
				}
				else if (dist1 < 4.82 && dist2 < 5.24){
					dihe=Vector::dihedral(atm1->getCoor(), atm2->getCoor(),atm3->getCoor(),atm4->getCoor());
					if (dihe > 42.1 && dihe < 119.5){
						atm1->setSS("G");
					}
				}
				else{
					//Do nothing
				}
			}
		}
	}

	for (i=0; i< ssmol->getAtmVecSize(); i++){
		std::cerr << Residue::aa321(ssmol->getAtom(i)->getResName()) << ssmol->getAtom(i)->getResId() << ":";
		std::cout << i+1 << ":";
		if (ssmol->getAtom(i)->getSS().length() == 1){
			std::cout << ssmol->getAtom(i)->getSS() << std::endl;
		}
		else{
			std::cout << "-" << std::endl;
		}
	}

	return ssmol;
}
