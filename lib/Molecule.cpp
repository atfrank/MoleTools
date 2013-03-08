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
  if (ifile.length() == 0){
    std::cerr << "Error: PDB file \"" << ifile << "\" cannot be found" << std::endl;
    return new Molecule;
  }
  else{
    return PDB::readPDB(ifile, model);
  }
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

//Analysis Functions

Vector Molecule::centerOfGeometry(bool selFlag){
	Vector cog=Vector(0.0, 0.0, 0.0);

	for (unsigned int i=0; i< this->getAtmVecSize(); i++){
		if (selFlag == true && this->getAtom(i)->getSel() == false){
			continue;
		}
		cog+=this->getAtom(i)->getCoor();
	}

	cog/=this->getAtmVecSize();

	return cog;
}

double Molecule::lsqfit (Molecule *refmol, bool move){
	Eigen::MatrixXd cmp; //Nx3 Mobile Coordinate Matrix
	Eigen::MatrixXd ref; //Nx3 Stationary Coordinate Matrix
	Eigen::Matrix3d R; //3x3 Covariance Matrix
	Eigen::Matrix3d C; //3x3 Covariance Matrix
	unsigned int i, j;
	Vector cmpCOG;
	Vector refCOG;
	Atom *atm;
	double E0;
	double RMSD;

	E0=0.0;
	
	//Need to resize dynamic matrix!
	cmp.resize(this->getAtmVecSize(),3);
	ref.resize(refmol->getAtmVecSize(),3);

	//Check selection sizes and load matrices

	cmpCOG=this->centerOfGeometry();

  for (i=0; i< this->getAtmVecSize(); i++){
    atm=this->getAtom(i);
    cmp.row(i) << atm->getX()-cmpCOG.x(), atm->getY()-cmpCOG.y(), atm->getZ()-cmpCOG.z();
		for (j=0; j< 3; j++){
			E0+=cmp(i,j)*cmp(i,j);
		}
  }

	refCOG=refmol->centerOfGeometry();

	for (i=0; i< refmol->getAtmVecSize(); i++){
	  atm=refmol->getAtom(i);
		ref.row(i) << atm->getX()-refCOG.x(), atm->getY()-refCOG.y(), atm->getZ()-refCOG.z();
		for (j=0; j< 3; j++){
      E0+=ref(i,j)*ref(i,j);
    }
	}

  //R=cmp.transpose()*ref;
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

 	if (move == true){
		Eigen::MatrixXd tmp;
		tmp=Umin*cmp.transpose();
		cmp=tmp.transpose();

		for (i=0; i< this->getAtmVecSize(); i++){
	  	atm=this->getAtom(i);
			atm->setCoor(Vector(cmp(i,0), cmp(i,1), cmp(i,2))+refCOG);
		}
  }

	RMSD=sqrt((E0-2.0*(S(0,0)+S(1,0)+S(2,0)))/this->getAtmVecSize());

	return RMSD;
}


