//Sean M. Law

#include "Molecule.hpp"

//Circular dependencies must be listed here instead of the header file
#include "PDB.hpp"
#include "Select.hpp"
#include "Analyze.hpp"

Molecule::~Molecule (){
	Chain *c;
	Residue *r;
	Atom *a;

	for (unsigned int i=0; i< this->getChnVecSize(); i++){
		c=this->getChain(i);
		for (unsigned int j=0; j< c->getResVecSize(); j++){
			r=c->getResidue(j);
			if (this->getCopyFlag() == false){
				//Only destruct atoms if molecule is NOT a copy
				for (unsigned int k=0; k< r->getAtmVecSize(); k++){
					a=r->getAtom(k);
					delete a;
				}
			}
			delete r;
		}
		delete c;
	}
}

Molecule* Molecule::readPDB (std::string ifile, int model){
  if (ifile.length() == 0){
    std::cerr << "Error: PDB file \"" << ifile << "\" cannot be found" << std::endl;
    return new Molecule;
  }
  else{
    return PDB::readPDB(ifile, model);
  }
}

std::string Molecule::writePDB(bool selFlag, bool print){

  std::ostringstream out;

  PDB::writePDBFormat(this, out, selFlag);

	if (print == true){
  	std::cout << out.str();
	}
	
	return out.str();
}

void Molecule::addAtom(Atom* atmEntry) {
  if (atmEntry->getAtmNum()){
    atmVec.push_back(atmEntry);
  }
}

Molecule* Molecule::clone (bool selFlag, bool keep){
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
	mol->setCopyFlag(false); //Not a copy
	
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

  if (keep == false){
    delete this;
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
	mol->setCopyFlag(true); //Is a copy, do not destruct atoms!

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
  for (unsigned int i=0; i< this->getNAtom(); i++){
    this->getAtom(i)->setSel(true);
  }
}

void Molecule::deselAll(){
  for (unsigned int i=0; i< this->getNAtom(); i++){
    this->getAtom(i)->setSel(false);
  }
}

void Molecule::select(std::string sel){
  Select::makeSel(this, sel);
}

unsigned int Molecule::getNAtom(){
  return atmVec.size();
}

unsigned int Molecule::getNAtomSelected(){
  unsigned int natom=0;
  for (unsigned int i=0; i< this->getAtmVecSize(); i++){
    if (this->getAtom(i)->getSel() == true){
      natom++;
    }
  }
  return natom;
}

void Molecule::setCopyFlag(bool copyFlagIn){
	copyFlag=copyFlagIn;
}

bool Molecule::getCopyFlag(){
	return copyFlag;
}

double Molecule::lsqfit (Molecule *refmol, bool transform){
	Eigen::MatrixXd cmp; //Nx3 Mobile Coordinate Matrix
	Eigen::MatrixXd ref; //Nx3 Stationary Coordinate Matrix
	Eigen::Matrix3d R; //3x3 Covariance Matrix
	unsigned int i;
	Vector cmpCOG;
	Vector refCOG;
	Vector coor;
	Atom *atm;
	//double E0;
	//unsigned int j;
	//double RMSD;
	unsigned int nrow;
	bool selFlag;

	//For optimal efficiency, atoms need to be pre-selected 
	//for both *this and *refmol before lsqfit() is called
	selFlag=true;

	//E0=0.0;		

	//Check selection sizes and resize matrices 
	if (this->getNAtomSelected() != refmol->getNAtomSelected()){
		std::cerr << "Error: Atom number mismatch in least squares fitting" << std::endl;
		return -1.0;
	}
	else{
		cmp.resize(this->getNAtomSelected(),3);
 		ref.resize(refmol->getNAtomSelected(),3);
	}

	cmpCOG=Analyze::centerOfGeometry(this, selFlag);

	nrow=0;
  for (i=0; i< this->getNAtom(); i++){
    atm=this->getAtom(i);
		if (atm->getSel() == false){
			continue;
		}
    cmp.row(nrow) << atm->getX()-cmpCOG.x(), atm->getY()-cmpCOG.y(), atm->getZ()-cmpCOG.z();
		//for (j=0; j< 3; j++){
		//	E0+=cmp(nrow,j)*cmp(nrow,j);
		//}
		nrow++;
  }

	refCOG=Analyze::centerOfGeometry(refmol, selFlag);

	nrow=0;
	for (i=0; i< refmol->getNAtom(); i++){
	  atm=refmol->getAtom(i);
		if (atm->getSel() == false){
			continue;
		}
		ref.row(nrow) << atm->getX()-refCOG.x(), atm->getY()-refCOG.y(), atm->getZ()-refCOG.z();
		//for (j=0; j< 3; j++){
    //  E0+=ref(nrow,j)*ref(nrow,j);
    //}
		nrow++;
	}

	R=ref.transpose()*cmp;

	Eigen::JacobiSVD<Eigen::MatrixXd> svd(R, Eigen::ComputeThinU | Eigen::ComputeThinV); 
	Eigen::MatrixXd V; 
	Eigen::MatrixXd S;
	Eigen::MatrixXd W;
	Eigen::MatrixXd Wt; //W transpose
	Eigen::MatrixXd Umin; //Rotation matrix
	double VdWtd; //Vt determinant * W determinant

	V=svd.matrixU();
	S=svd.singularValues();
	W=svd.matrixV();

	Wt=W.transpose();

	VdWtd=V.determinant()*Wt.determinant();

  if (VdWtd < 0.0){
    //Is a reflection!
    //Make third column negative
    //Untested
    std::cerr << "Warning: LSQFIT Relection Found!" << std::endl;
    for (i=0; i< 3; i++){
      V(i,2)=-V(i,2); //Base zero
    }
    S(2,0)=-S(2,0);
  }

	Umin=V*Wt; //Optimal rotation matrix, should already be normalized

 	if (transform == true){
		//Apply best fit (rigid body) transformation to entire molecule
		for (i=0; i< this->getNAtom(); i++){
	  	atm=this->getAtom(i);
			atm->setCoor(atm->getCoor()-cmpCOG); //Translate molecule to origin based on selected
			//Apply rotation matrix to entire molecule, not just selected
			//This is equivalent to doing Umin*cmp but for the entire molecule, not just selected
			coor.x()=Umin(0,0)*atm->getX()+Umin(0,1)*atm->getY()+Umin(0,2)*atm->getZ();
			coor.y()=Umin(1,0)*atm->getX()+Umin(1,1)*atm->getY()+Umin(1,2)*atm->getZ();
			coor.z()=Umin(2,0)*atm->getX()+Umin(2,1)*atm->getY()+Umin(2,2)*atm->getZ();
			coor+=refCOG; //Translate molecule to reference based on selected
			atm->setCoor(coor);
		}
  }

	//RMSD=sqrt((E0-2.0*(S(0,0)+S(1,0)+S(2,0)))/this->getNAtom());
	//std::cerr << "RMSD of selected : " << RMSD << std::endl;

	return S(0,0)+S(1,0)+S(2,0);
}

double Molecule::rmsd (Molecule *refmol){
	return Analyze::rmsd(this, refmol);
}

void Molecule::recenter (Molecule *recmol){
	std::cerr << "Warning: Image recentering has not been implemented yet" << std::endl;
}

void Molecule::translate (const double &dx, const double &dy, const double &dz){
  for(unsigned int i=0; i< this->getNAtom(); i++){
    this->getAtom(i)->setCoor(this->getAtom(i)->getCoor()+Vector(dx, dy, dz));
  }
}

void Molecule::translate (const Vector &u){
  for(unsigned int i=0; i< this->getNAtom(); i++){
    this->getAtom(i)->setCoor(this->getAtom(i)->getCoor()+u);
  }
}

void Molecule::rotate (const double &r1c1, const double &r1c2, const double &r1c3, const double &r2c1, const double &r2c2, const double &r2c3, const double &r3c1, const double &r3c2, const double &r3c3){
  Atom *atm;
  Vector coor;
  Vector cog;

  cog=Analyze::centerOfGeometry(this, false);

  for (unsigned int i=0; i< this->getNAtom(); i++){
    atm=this->getAtom(i);
    atm->setCoor(atm->getCoor()-cog); //Translate molecule to origin based on selected
    //Apply rotation matrix to entire molecule
    coor.x()=r1c1*atm->getX()+r1c2*atm->getY()+r1c3*atm->getZ();
    coor.y()=r2c1*atm->getX()+r2c2*atm->getY()+r2c3*atm->getZ();
    coor.z()=r3c1*atm->getX()+r3c2*atm->getY()+r3c3*atm->getZ();
    coor+=cog; //Translate molecule back
    atm->setCoor(coor);
  }
}

void Molecule::center (bool selFlag){
  Vector cog=Analyze::centerOfGeometry(this, selFlag);
  for (unsigned int i=0; i< this->getNAtom(); i++){
    this->getAtom(i)->setCoor(this->getAtom(i)->getCoor()-cog);
  }
}
