//Sean M. Law

#include "Analyze.hpp"

Analyze::Analyze (){
	type.clear();
	sel.clear();
	mol.clear();
	tdata.clear();
	ndata=0;
	resel=false;
}

void Analyze::setType(const std::string& typein){
	this->type=typein;
	//Depending upon the type, initialize necessary vectors
}

std::string Analyze::getType(){
	return this->type;
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
	unsigned int i;
	Molecule* tmpmol;

	if (this->getType() == "newtype"){

	}
	else if (this->getType() == "rmsf"){
    molin->select(this->getSel(0));
    tmpmol=molin->copy();
    this->addMol(tmpmol);
  }
	else if (this->getType() == "rmsd"){
		molin->select(this->getSel(0));
		tmpmol=molin->copy();
		this->addMol(tmpmol);
	}
	else if (this->getType() == "average"){
		molin->select(this->getSel(0));
		tmpmol=molin->copy();
		this->addMol(tmpmol);
	}
	else{
		for (i=0; i< this->getNSel(); i++){
			molin->select(this->getSel(i));
			tmpmol=molin->copy();
			this->addMol(tmpmol);
		}
	}
}

Molecule* Analyze::getMol(const int& element){
	return this->mol.at(element);
}

unsigned int Analyze::getNMol(){
	return this->mol.size();
}

void Analyze::runAnalysis(){
	Vector xyz;

	std::cout << std::fixed;
	//Process all analyses based on type
	if (this->getType() == "newtype"){
		
	}
	else if (this->getType() == "rmsf"){
		//Resize and initialize if necessary
		if (ndata == 0){
			tdata.resize(this->getMol(0)->getNAtomSelected(), 0);
		}
		//No output during analysis. Only during post-analysis
		Analyze::rmsf(this->getMol(0), this->getMol(1), tdata, ndata);
	}
	else if (this->getType() == "rmsd"){
		if (this->getNMol() == 2){
			std::cout << std::setw(9) << std::right << std::setprecision(3) << Analyze::rmsd(this->getMol(0),this->getMol(1));
		}
		else{
			std::cout << std::setw(9) << std::right << std::setprecision(3) << "NaN";
		}
	}
	else if (this->getType() == "average"){
		//No output during analysis. Only during post-analysis
		Analyze::avgMol(this->getMol(0), this->getMol(1), ndata);
	}
	else if(this->getType() == "quick"){
		if (this->getNMol() == 1){
			xyz=Analyze::centerOfGeometry(this->getMol(0));
     	std::cout << std::setw(9) << std::right << std::setprecision(3) << xyz.x();
			std::cout << std::setw(9) << std::right << std::setprecision(3) << xyz.y();
			std::cout << std::setw(9) << std::right << std::setprecision(3) << xyz.z();
		}
		else if (this->getNMol() == 2){
      std::cout << std::setw(9) << std::right << std::setprecision(3) << Analyze::distance(Analyze::centerOfGeometry(mol.at(0)), Analyze::centerOfGeometry(mol.at(1)));
		}
		else if (this->getNMol() == 3){
			std::cout << std::setw(9) << std::right << std::setprecision(3) << Analyze::angle(Analyze::centerOfGeometry(this->getMol(0)), Analyze::centerOfGeometry(this->getMol(1)), Analyze::centerOfGeometry(this->getMol(2)));
		}
		else if (this->getNMol() == 4){
			std::cout << std::setw(9) << std::right << std::setprecision(3) << Analyze::dihedral(Analyze::centerOfGeometry(this->getMol(0)), Analyze::centerOfGeometry(this->getMol(1)), Analyze::centerOfGeometry(this->getMol(2)), Analyze::centerOfGeometry(this->getMol(3)));
		}
		else{
			std::cout << std::fixed;
      std::cout << std::setw(9) << "NaN";
		}
	}
	else{
		std::cerr << "Error: Skipping unrecognized analysis type \"" << type << "\"" << std::endl;
	}
}

void Analyze::postAnalysis(){
	unsigned int i, j;
	Atom *atm;

	std::cout << std::fixed;
	std::cout << std::endl << std::endl; //For using index in Gnuplot
  //Process all analyses based on type
  if (this->getType() == "rmsf"){
		j=0;
		for(i=0; i< this->getMol(0)->getNAtomSelected(); i++){
			atm=this->getMol(0)->getAtom(i);
			if (atm->getSel() == false){
				continue;
			}
			tdata.at(j)=sqrt(tdata.at(j)/ndata);
			std::cout << atm->getSummary();  
			std::cout << std::setw(9) << std::right << std::setprecision(3) << tdata.at(j);
			std::cout << std::endl;
			j++;
		}
	}
	else if (this->getType() == "average"){
		for(i=0; i< this->getMol(1)->getNAtomSelected(); i++){
			atm=this->getMol(1)->getAtom(i);
      if (atm->getSel() == false){
        continue;
      }
			atm->setCoor(atm->getCoor()/ndata);
		}
		this->getMol(1)->writePDB(true, true);
	}
	else{
		//Do nothing
	}
}

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

void Analyze::avgMol (Molecule* cmpmol, Molecule* refmol, int &ndataIO){
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
