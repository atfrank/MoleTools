//Sean M. Law

#include "Analyze.hpp"

Analyze::Analyze (){
	sel.clear();
	mol.clear();
	tdata.clear();
	ndata=0;
	resel=false;
}

void Analyze::addSel(const std::string& selin){
	this->sel.push_back(selin);
}

std::string Analyze::getSel(const int& element){
	return this->sel.at(element);
}

unsigned int Analyze::getNSel(){
	return this->sel.size();
}

void Analyze::addMol(Molecule* molin){
	this->mol.push_back(molin);
}

void Analyze::setMol(const int element, Molecule* molin){
	this->mol.at(element)=molin;
}

void Analyze::resizeNMol(const int sizein){
  this->mol.resize(sizein);
}

void Analyze::setupMolSel(Molecule* molin){
	Molecule* tmpmol;

  molin->select(this->getSel(0));
  tmpmol=molin->copy();
  this->addMol(tmpmol);
}

void AnalyzeDistance::setupMolSel(Molecule* molin){
	unsigned int i;
  Molecule* tmpmol;

	for (i=0; i< this->getNSel(); i++){
		molin->select(this->getSel(i));
		tmpmol=molin->copy();
		this->addMol(tmpmol);
	}
}

void AnalyzeAngle::setupMolSel(Molecule* molin){
  unsigned int i;
  Molecule* tmpmol;

  for (i=0; i< this->getNSel(); i++){
    molin->select(this->getSel(i));
    tmpmol=molin->copy();
    this->addMol(tmpmol);
  }
}

void AnalyzeDihedral::setupMolSel(Molecule* molin){
  unsigned int i;
  Molecule* tmpmol;

  for (i=0; i< this->getNSel(); i++){
    molin->select(this->getSel(i));
    tmpmol=molin->copy();
    this->addMol(tmpmol);
  }
}

Molecule* Analyze::getMol(const int& element){
	return this->mol.at(element);
}

unsigned int Analyze::getNMol(){
	return this->mol.size();
}

void Analyze::setNData(const int& ndatain){
	ndata=ndatain;
}

int& Analyze::getNData(){
	return ndata;
}

std::vector<double>& Analyze::getTDataVec(){
	return tdata;
}

//All preAnalysis Functions

void Analyze::preAnalysis(Molecule* molin){
	this->setupMolSel(molin);
}

void AnalyzeRMSD::preAnalysis(Molecule* molin){
	this->setupMolSel(molin);
	Molecule* refmol=molin->clone();
	this->setupMolSel(refmol);
}

void AnalyzeRMSF::preAnalysis(Molecule* molin){
  this->setupMolSel(molin);
  Molecule* refmol=molin->clone();
  this->setupMolSel(refmol);
}

void AnalyzeAverage::preAnalysis(Molecule* molin){
  this->setupMolSel(molin);
  Molecule* refmol=molin->clone();
  this->setupMolSel(refmol);
}


//All runAnalysis Functions

void AnalyzeCOG::runAnalysis(){
  Vector xyz=Analyze::centerOfGeometry(this->getMol(0));
  std::cout << std::fixed;
  std::cout << std::setw(9) << std::right << std::setprecision(3) << xyz.x();
  std::cout << std::setw(9) << std::right << std::setprecision(3) << xyz.y();
  std::cout << std::setw(9) << std::right << std::setprecision(3) << xyz.z();
}

void AnalyzeRMSD::runAnalysis(){
  std::cout << std::fixed;
  if (this->getNMol() == 2){
    std::cout << std::setw(9) << std::right << std::setprecision(3) << Analyze::rmsd(this->getMol(0),this->getMol(1));
  }
  else{
    std::cout << std::setw(9) << std::right << std::setprecision(3) << "NaN";
  }
}

void AnalyzeRMSF::runAnalysis(){
  std::vector<double> &tdata=this->getTDataVec();

  //Resize and initialize if necessary
  if (getNData() == 0){
    tdata.resize(this->getMol(0)->getNAtomSelected(), 0);
  }
  Analyze::rmsf(this->getMol(0), this->getMol(1), tdata, getNData());
}

void AnalyzeAverage::runAnalysis(){
  if (getNData() == 0){
    this->getMol(1)->zeroCoor();
  }
  Analyze::averageMol(this->getMol(0), this->getMol(1), getNData());
}

void AnalyzeDistance::runAnalysis(){
	std::cout << std::fixed;
	std::cout << std::setw(9) << std::right << std::setprecision(3) << Analyze::distance(Analyze::centerOfGeometry(this->getMol(0)), Analyze::centerOfGeometry(this->getMol(1)));
}

void AnalyzeAngle::runAnalysis(){
	std::cout << std::fixed;
  std::cout << std::setw(9) << std::right << std::setprecision(3) << Analyze::angle(Analyze::centerOfGeometry(this->getMol(0)), Analyze::centerOfGeometry(this->getMol(1)), Analyze::centerOfGeometry(this->getMol(2)));
}

void AnalyzeDihedral::runAnalysis(){
	std::cout << std::fixed;
  std::cout << std::setw(9) << std::right << std::setprecision(3) << Analyze::dihedral(Analyze::centerOfGeometry(this->getMol(0)), Analyze::centerOfGeometry(this->getMol(1)), Analyze::centerOfGeometry(this->getMol(2)), Analyze::centerOfGeometry(this->getMol(3)));
}

//All postAnalysis functions

void Analyze::postAnalysis(){
	//Do nothing
}

void AnalyzeAverage::postAnalysis(){
  unsigned int i;
  Atom *atm;

  std::cout << std::fixed;
  std::cout << std::endl << std::endl; //For using index in Gnuplot
  for(i=0; i< this->getMol(1)->getNAtomSelected(); i++){
    atm=this->getMol(1)->getAtom(i);
    if (atm->getSel() == false){
      continue;
    }
    atm->setCoor(atm->getCoor()/getNData());
  }
  this->getMol(1)->writePDB();
}

void AnalyzeRMSF::postAnalysis(){
  unsigned int i, j;
  Atom *atm;
  std::vector<double> &tdata=this->getTDataVec();

  std::cout << std::fixed;
  std::cout << std::endl << std::endl; //For using index in Gnuplot
  j=0;
  for(i=0; i< this->getMol(0)->getNAtomSelected(); i++){
    atm=this->getMol(0)->getAtom(i);
    if (atm->getSel() == false){
      continue;
    }
    tdata.at(j)=sqrt(tdata.at(j)/getNData());
    std::cout << atm->getSummary();
    std::cout << std::setw(9) << std::right << std::setprecision(3) << tdata.at(j);
    std::cout << std::endl;
    j++;
  }
}



//Basic analysis functions

Vector Analyze::centerOfGeometry(Molecule *mol, bool selFlag){
	Vector cog=Vector(0.0, 0.0, 0.0);

  for (unsigned int i=0; i< mol->getAtmVecSize(); i++){
    if (selFlag == true && mol->getAtom(i)->getSel() == false){
      continue;
    }
    cog+=mol->getAtom(i)->getCoor();
  }

  if (selFlag == true){
    cog/=mol->getNAtomSelected();
  }
  else{
    cog/=mol->getNAtom();
  }

	return cog;
}

double Analyze::rmsd (Molecule *cmpmol, Molecule *refmol){
  double RMSD;
  double E;
  unsigned int i, j;
  Atom *atm;
  std::vector<double> dx;
  std::vector<double> dy;
  std::vector<double> dz;

  //For optimal efficiency, atoms need to be pre-selected
  //for both *cmpmol and *refmol before rmsd() is called
  RMSD=-1.0;
  E=0.0;

  //Check selection sizes and resize matrices
  if (cmpmol->getNAtomSelected() != refmol->getNAtomSelected()){
    std::cerr << std::endl << "Error: Atom number mismatch in RMSD calculation" << std::endl;
		std::cerr << "CMP-NATOM: " << cmpmol->getNAtomSelected() << ", ";
		std::cerr << "REF-NATOM: " << refmol->getNAtomSelected() << std::endl;
    return -1.0;
  }

  for (i=0; i< cmpmol->getNAtom(); i++){
    atm=cmpmol->getAtom(i);
    if (atm->getSel() == false){
      continue;
    }
    dx.push_back(atm->getX());
    dy.push_back(atm->getY());
    dz.push_back(atm->getZ());
  }

  j=0;
  for (i=0; i< refmol->getNAtom(); i++){
    atm=refmol->getAtom(i);
    if (atm->getSel() == false){
      continue;
    }
    dx.at(j)=dx.at(j)-atm->getX();
    dy.at(j)=dy.at(j)-atm->getY();
    dz.at(j)=dz.at(j)-atm->getZ();
    E+=dx.at(j)*dx.at(j)+dy.at(j)*dy.at(j)+dz.at(j)*dz.at(j);
    j++;
  }

  RMSD=sqrt(E/cmpmol->getNAtomSelected());

  return RMSD;

}

void Analyze::rmsf (Molecule* cmpmol, Molecule* refmol, std::vector<double> &tdataIO, int &ndataIO){

	unsigned i,j;
	Atom *atm;
	std::vector<Vector> coor;
	Vector d;

	//Check selection sizes and resize matrices
  if (cmpmol->getNAtomSelected() != refmol->getNAtomSelected()){
		std::cerr << std::endl << "Error: Atom number mismatch in RMSD calculation" << std::endl;
    std::cerr << "CMP-NATOM: " << cmpmol->getNAtomSelected() << ", ";
    std::cerr << "REF-NATOM: " << refmol->getNAtomSelected() << std::endl;
	}
	else{
		for (i=0; i< cmpmol->getNAtom(); i++){
			atm=cmpmol->getAtom(i);
    	if (atm->getSel() == false){
      	continue;
    	}
			coor.push_back(atm->getCoor());
		}
		
		j=0;
		for (i=0; i< refmol->getNAtom(); i++){
			atm=refmol->getAtom(i);
	    if (atm->getSel() == false){
	      continue;
	    }
			d=coor.at(j)-atm->getCoor();
			tdataIO.at(j)+=d.norm()*d.norm(); //Add squared distance to previous frame
			j++;
		}
		ndataIO++;
	}
}

void Analyze::averageMol (Molecule* cmpmol, Molecule* refmol, int &ndataIO){
  unsigned i,j;
  Atom *atm;
  std::vector<Vector> coor;
  Vector d;

  //Check selection sizes and resize matrices
  if (cmpmol->getNAtomSelected() != refmol->getNAtomSelected()){
    std::cerr << std::endl << "Error: Atom number mismatch in RMSD calculation" << std::endl;
    std::cerr << "CMP-NATOM: " << cmpmol->getNAtomSelected() << ", ";
    std::cerr << "REF-NATOM: " << refmol->getNAtomSelected() << std::endl;
  }
  else{
    for (i=0; i< cmpmol->getNAtom(); i++){
      atm=cmpmol->getAtom(i);
      if (atm->getSel() == false){
        continue;
      }
      coor.push_back(atm->getCoor());
    }

    j=0;
    for (i=0; i< refmol->getNAtom(); i++){
      atm=refmol->getAtom(i);
      if (atm->getSel() == false){
        continue;
      }
      atm->setCoor(coor.at(j)+atm->getCoor());
      j++;
    }
    ndataIO++;
  }
}

double Analyze::distance (const Vector& u, const Vector& v){
	Vector d=u-v;
	return d.norm();
}

double Analyze::angle (const Vector& u, const Vector& v, const Vector& w){
  double angle;
  Vector dx, dy;
  double dp, nx, ny;

  dx=u-v;
  dy=w-v;

  nx=dx.norm();
  ny=dy.norm();

  dp=dx.dot(dy);

  angle=acos(dp/(nx*ny));

  return angle/PI*180.0;
}

double Analyze::dihedral (const Vector& t, const Vector& u, const Vector& v, const Vector& w) {
  double dihedral;
  Vector dx, dy, dz, p1, p2, p3;
  double np1, np2, dp1, dp2, ts;

  dx=t-u;
  dy=u-v;
  dz=w-v; //This is correct!

  p1=dx.cross(dy);

  np1=p1.norm();
  p1=p1/np1;

  p2=dz.cross(dy);
  np2=p2.norm();
  p2=p2/np2;

  dp1=p1.dot(p2); //Dot product

  ts=1.0-dp1*dp1;
  ts=(ts<0.0)?0.0:sqrt(ts);
  dihedral=PI/2.0-atan2(dp1,ts);

  p3=p1.cross(p2);

  dp2=p3.dot(dy); //Dot product

  if (dp2 > 0.0){
    dihedral=-dihedral;
  }

  return dihedral/PI*180.0;
}

double Analyze::distance (Molecule* sel1, Molecule* sel2, bool selFlag){
	return Analyze::distance(Analyze::centerOfGeometry(sel1,selFlag), Analyze::centerOfGeometry(sel2,selFlag));
}

double Analyze::angle (Molecule* sel1, Molecule* sel2, Molecule* sel3, bool selFlag){
	return Analyze::angle(Analyze::centerOfGeometry(sel1,selFlag), Analyze::centerOfGeometry(sel2,selFlag), Analyze::centerOfGeometry(sel3,selFlag));
}

double Analyze::dihedral (Molecule* sel1, Molecule* sel2, Molecule* sel3, Molecule* sel4, bool selFlag){
	return Analyze::dihedral(Analyze::centerOfGeometry(sel1,selFlag), Analyze::centerOfGeometry(sel2,selFlag), Analyze::centerOfGeometry(sel3,selFlag), Analyze::centerOfGeometry(sel4,selFlag));
}

void Analyze::pairwiseDistance(Molecule *mol, std::map<std::pair<Atom*, Atom*>, double>& pdin){
	unsigned int i, j;
	double distance;

	pdin.clear();

	for (i=0; i< mol->getAtmVecSize(); i++){
		pdin.insert(std::make_pair(std::make_pair(mol->getAtom(i), mol->getAtom(i)), 0.0)); //Zero diagonal
	}

  for (i=0; i< mol->getAtmVecSize(); i++){
    for (j=i+1; j< mol->getAtmVecSize(); j++){
			distance=9999.9;
			if (mol->getAtom(i)->getX() < 9999.9 && mol->getAtom(j)->getX() < 9999.9){
				distance=Analyze::distance(mol->getAtom(i)->getCoor(), mol->getAtom(j)->getCoor());
			}
			 pdin.insert(std::make_pair(std::make_pair(mol->getAtom(i), mol->getAtom(j)), distance));
			  pdin.insert(std::make_pair(std::make_pair(mol->getAtom(j), mol->getAtom(i)), distance));
    }
  }
}

void Analyze::allAnglesDihedrals(Molecule *mol, std::map<Atom*, std::pair<double, double> >& adin){
	unsigned int i, j;
	Chain *c;
	Atom *atm1, *atm2, *atm3, *atm4;
	unsigned int size;
	double angle, dihedral;

	for (i=0; i< mol->getChnVecSize(); i++){
    c=mol->getChain(i);
    size=c->getAtmVecSize();
    for (j=0; j< c->getAtmVecSize(); j++){
      atm1=c->getAtom(j);
      atm2=NULL;
      atm3=NULL;
      atm4=NULL;

			if (j+1 < size && atm1->getResId()+1 == c->getAtom(j+1)->getResId()){
				atm2=c->getAtom(j+1);
			}
			else{
				if (j+1 < size && atm1->getResId() == c->getAtom(j+1)->getResId() && (atm1->getICode().compare(0,1,c->getAtom(j+1)->getICode(),0,1) != 0)){
					atm2=c->getAtom(j+1);
				}
			}

			if (j+2 < size && atm1->getResId()+2 == c->getAtom(j+2)->getResId()){
				atm3=c->getAtom(j+2);
			}
		 	else{
        if (j+2 < size && atm1->getResId() == c->getAtom(j+2)->getResId() && (atm1->getICode().compare(0,1,c->getAtom(j+2)->getICode(),0,1) != 0)){
          atm3=c->getAtom(j+2);
        }
      }

			if (j+3 < size && atm1->getResId()+3 == c->getAtom(j+3)->getResId()){
        atm4=c->getAtom(j+3);
      }
			else{
        if (j+3 < size && atm1->getResId() == c->getAtom(j+3)->getResId() && (atm1->getICode().compare(0,1,c->getAtom(j+3)->getICode(),0,1) !=0)){
          atm4=c->getAtom(j+3);
        }
      }

			angle=9999.9;
			dihedral=9999.9;
			if (atm2 != NULL && atm3 != NULL){
				//Get Angle
				angle=Analyze::angle(atm1->getCoor(), atm2->getCoor(), atm3->getCoor());
				if (atm4 != NULL){
					//Get Dihedral
					dihedral=Analyze::dihedral(atm1->getCoor(), atm2->getCoor(), atm3->getCoor(), atm4->getCoor());
				}
			}
			adin.insert(std::make_pair(atm1, std::make_pair(angle, dihedral)));
    }
  }
}

void Analyze::pcasso(Molecule* mol){
	Chain *c;
	Atom *ai, *aj;
	unsigned int j, start;
	double defVal;
	unsigned int icapc; //Counter for Ca or pseudocenter (pc)
	double iMinus6, iPlus6; //For shortest non-local contact
	std::map<std::pair<Atom*, Atom*>, double> caPairDist; //Ca-Ca Distances
	std::map<std::pair<Atom*, Atom*>, double> pcPairDist; //Pseudocenter-Pseudocenter Distances
	std::map<Atom*, std::pair<double, double> > caAngDihe; //Ca-Ca Angle/Diehdral
	std::map<Atom*, std::pair<double, double> > pcAngDihe; //Pseudocenter-Pseudocenter Angle/Dihedral
	Molecule* cmol;

	defVal=9999.9;

	mol->storeSel();
	mol->select(":.CA");
	cmol=mol->clone(true,true); //Copy selection, keep original
	mol->recallSel(); //Restore original selection
	mol->eraseSel();

	//Analyze all C-alpha first
	Analyze::pairwiseDistance(cmol, caPairDist);
	Analyze::allAnglesDihedrals(cmol, caAngDihe);

	//Analyze all pseudocenter
	cmol->modPseudoCenter();
	Analyze::pairwiseDistance(cmol, pcPairDist);
	Analyze::allAnglesDihedrals(cmol, pcAngDihe);

	for (unsigned int ichain=0; ichain < cmol->getChnVecSize(); ichain++){
		c=cmol->getChain(ichain);
		for (unsigned int iatom=0; iatom < c->getAtmVecSize(); iatom++){
			ai=c->getAtom(iatom);
			std::cout << ai->getPdbId() << " " << ai->getChainId() << " " << ai->getResId() << " ";
			//i-2, i-1, i+1, i+2, i+3, i+4, i+5 Distances
			for (icapc=0; icapc<=1; icapc++){
				//Deal with unsigned int subtraction from zero
				if (iatom == 0){
					std::cout << defVal << " " << defVal << " ";
					start=iatom-0;
				}
				else if (iatom == 1){
					std::cout << defVal << " ";
					start=iatom-1;
				}
				else{
					start=iatom-2;
				}
				for (j=start; j<= iatom+5; j++){
					aj=c->getAtom(j);
					if (aj == NULL){
						std::cout << defVal << " ";
					}
					else if (j == iatom){
						//Distance == 0
						continue;
					}
					else{
						if (icapc == 0){ //Ca
							std::cout << caPairDist.at(std::make_pair(ai, aj)) << " ";
						}
						else{ //Pseudocenter
							std::cout << pcPairDist.at(std::make_pair(ai, aj)) << " ";
						}
					}
				}
			}

			//Angles and Dihedrals
			for (icapc=0; icapc<=1; icapc++){
				if (icapc == 0){ //Ca
					std::cout << caAngDihe.at(ai).first << " " << caAngDihe.at(ai).second << " ";
				}
				else{ //Pseudocenter
					std::cout << pcAngDihe.at(ai).first << " " << pcAngDihe.at(ai).second << " ";
				}
			}

			//Shortest non-local contact distance, > i+6 and < i-6
			iPlus6=9999.9;
			iMinus6=9999.9;
			for (j=0; j< cmol->getAtmVecSize(); j++){
				aj=cmol->getAtom(j);
				if (ai == aj){
					continue;
				}
				if (caPairDist.at(std::make_pair(ai, aj)) < iPlus6){
					if (ai->getChainId().compare(aj->getChainId()) != 0)){
						//atom i and atom j are on different chains
						iPlus6=caPairDist.at(std::make_pair(ai, aj));
					}
					else{
						//atom i and atom j are on the same chain
						
					}
				}
				if (caPairDist.at(std::make_pair(ai, aj)) < iMinus6){
					if (ai->getChainId().compare(aj->getChainId()) != 0){
						//atom i and atom j are on different chains
						iMinus6=caPairDist.at(std::make_pair(ai, aj));
					}
					else{
						//atom i and atom j are on the same chain

					}
				}
			}
			std::cout << iPlus6 << " " << iMinus6 << " ";
			std::cout << std::endl;
		}
	}
}
