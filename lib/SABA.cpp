//Sean M. Law

#include "SABA.hpp"

Molecule* SABA::getPseudoCenter(Molecule *mol){
	Molecule *ssmol=new Molecule;
	Atom *atmEntry;
	Atom *lastAtom;
	Vector xyz;
	unsigned int i, j;
  std::string sel;

	//Leave original molecule untouched
  sel=mol->getChain(0)->getChainId();
  sel+=":.CA";
	ssmol=mol->clone();
	ssmol->select(sel);
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
	double dist1, dist2, dihe;

	//Alpha and 310 Helices
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
            //310 Helix
						atm1->setSS("G");
					}
				}
				else{
					//Do nothing
				}
			}
		}
	}

  size=ssmol->getAtmVecSize();
  for (i=1; i< size; i++){ //i can't be zero, no i-1
    atm1=ssmol->getAtom(i);
    atm2=ssmol->getAtom(i-1);
    if (atm1->getChainId() == atm2->getChainId() && atm1->getResId()-1 == atm2->getResId()){
      //Got i' and i'-1
      for (j=1; j< size-4; j++){ //j can't be zero, no j-1
        atm3=ssmol->getAtom(j);
        if (atm1->getChainId() == atm3->getChainId() && atm1->getResId() < atm3->getResId()+4 && atm1->getResId() > atm3->getResId()-4){
          continue;
        }
        atm4=ssmol->getAtom(j-1);
        if (atm3->getChainId() == atm4->getChainId() && atm3->getResId()-1 == atm4->getResId()){
          //Got j' and j'-1
          //Calculate i' and j' distance
          dist1=Vector::distance(atm1->getCoor(), atm3->getCoor());
          if (dist1 > 2.58 && dist1 <5.18){
            dist2=Vector::distance(atm2->getCoor(),atm4->getCoor());
            if (dist2 > 4.34 && dist2 < 5.03){
              atm1->setSS("E");
              atm2->setSS("E");
              atm3->setSS("E");
              atm4->setSS("E");
            }
          }
        }
      }
    }
    else if (i < size-1){
      //Anti-parallel
      //atm1 = i above
      atm2=ssmol->getAtom(i+1);
      if (atm1->getChainId() == atm2->getChainId() && atm1->getResId()+1 == atm2->getResId()){
        //Got i' and i'+1
        for (j=1; j< size-4; j++){ //j can't be zero, no j-1
          atm3=ssmol->getAtom(j);
          if (atm1->getChainId() == atm3->getChainId() && atm1->getResId() < atm3->getResId()+4 && atm1->getResId() > atm3->getResId()-4){
            continue;
          }
          atm4=ssmol->getAtom(j-1);
          if (atm3->getChainId() == atm4->getChainId() && atm3->getResId()-1 == atm4->getResId()){
            //Got j' and j'-1
            dist1=Vector::distance(atm1->getCoor(), atm3->getCoor());
            if (dist1 > 4.36 && dist1 < 5.19){
              dist2=Vector::distance(atm2->getCoor(), atm4->getCoor());
              if (dist2 > 4.16 && dist2 < 5.27){
                atm1=mol->getAtom(i+1);//C-alpha
                atm3=mol->getAtom(j); //C-alpha
                dist1=Vector::distance(atm1->getCoor(), atm2->getCoor());
                if (dist1 > 1.42 && dist1 < 5.99){
                  //Check exceptions
                  atm1->setSS("E");
                }
              }
            }
          }
        }
      }
    }
  }

	for (i=0; i< ssmol->getAtmVecSize(); i++){
		std::cout << ssmol->getAtom(i)->getChainId() << ":" <<  Residue::aa321(ssmol->getAtom(i)->getResName()) << ssmol->getAtom(i)->getResId(); 
		std::cout << "(" << i+1 << ")";
		if (ssmol->getAtom(i)->getSS().length() == 1){
			std::cout << ssmol->getAtom(i)->getSS() << std::endl;
		}
		else{
			std::cout << "-" << std::endl;
		}
	}

	return ssmol;
}
