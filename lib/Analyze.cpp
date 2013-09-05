//Sean M. Law

#include "Analyze.hpp"

Analyze::Analyze (){
	sel.clear();
	mol.clear();
	tdata.clear();
	ndata=0;
	resel=false;
  avgCovar.resize(0,0);
  ifile.clear();
  ofile.clear();
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

Eigen::MatrixXd& Analyze::getAvgCovar(){
  return avgCovar;
}

void Analyze::addModes(const std::vector<unsigned int>& modesin){
  modes=modesin;
}

void Analyze::initCovar(const unsigned int &xin, const unsigned int &yin){
  avgCovar=Eigen::MatrixXd::Zero(xin, yin);
}

void Analyze::diagonalizeCovar(){
  //Diagonalize avg covariance matrix
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigOut(avgCovar);
  //Eigenvalues are already sorted in increasing order, largest at N
  //std::cout << eigOut.eigenvalues() << std::endl;

  //Eigenvectors are stored in columns, starting with
  //(x1, y1, z1, x2, y2, z2, ..., xN, yN, zN)
  //std::cout << eigOut.eigenvectors() << std::endl;

  //Output PC1
  //unsigned int pc1 = 3*this->getMol(1)->getNAtomSelected() - 1;
  //std::cout << eigOut.eigenvalues()[pc1] << std::endl;
  //std::cout << eigOut.eigenvectors().col(pc1) << std::endl;

  //Output PC2
  //unsigned int pc2 = 3*this->getMol(1)->getNAtomSelected() - 2;
  //std::cout << eigOut.eigenvalues()[pc2] << std::endl;
  //std::cout << eigOut.eigenvectors().col(pc2) << std::endl;
}

void Analyze::setInput(const std::string& fin){
  ifile=fin;
}

std::string Analyze::getInput(){
  return ifile;
}

void Analyze::setOutput (const std::string& fin){
  ofile=fin;
}
std::string Analyze::getOutput(){
  return ofile;
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
  //Zero coordinates for future appending
  this->getMol(1)->zeroCoor();
}

void AnalyzeEDA::preAnalysis(Molecule* molin){
  this->setupMolSel(molin);
  Molecule* refmol=molin->clone();
  this->setupMolSel(refmol);
  //Resize matrix to 3N x 3N and zero
  this->initCovar(3*refmol->getNAtomSelected(), 3*refmol->getNAtomSelected());
}

void AnalyzeProjection::preAnalysis(Molecule* molin){
  std::ifstream inpFile;
  std::istream* inp;
  std::string line;
  unsigned int nrow;
  std::vector<double> tmp;
  unsigned int N3;

  this->setupMolSel(molin);
  Molecule* refmol=molin->clone();
  this->setupMolSel(refmol);

  //Resize matrix to 3N x 3N and zero
  N3=3*refmol->getNAtomSelected();
  this->initCovar(N3, N3);

  //Read covariance matrix
  if (getInput().length() > 0){
    nrow=0;

    inpFile.open(getInput().c_str(), std::ios::in);
    inp=&inpFile;

    while (inp->good() && !(inp->eof())){
      getline(*inp, line);
      if (line.length() > 0){
        Misc::splitNum(Misc::trim(line), " \t", tmp, false);
        if (tmp.size() != N3){
          std::cerr << "Warning: The number of columns in the covariance matrix (" << tmp.size();
          std::cerr << ") does not match 3x the number of selected atoms (" << N3;
          std::cerr << ")" << std::endl;
        }
        for (unsigned int ncol=0; ncol< N3 && ncol < tmp.size(); ncol++){
          getAvgCovar()(nrow, ncol)=tmp.at(ncol);
        }
        nrow++;
      }
    }
    if (nrow != N3){
      std::cerr << "Warning: The number of rows in the covariance matrix (" << nrow;
      std::cerr << ") does not match 3x the number of selected atoms (" << N3;
      std::cerr << ")" << std::endl;
    }

    if (inpFile.is_open()){
      inpFile.close();
    }
  }

  std::cout << getAvgCovar() << std::endl;
  //Diagonalize covariance matrix
  diagonalizeCovar();
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

void AnalyzeEDA::runAnalysis(){
  Analyze::averageCovariance(this->getMol(0), this->getMol(1), getAvgCovar(), getNData());
}

void AnalyzeProjection::runAnalysis(){
  Analyze::projectModes(this->getMol(0), this->getMol(1), getAvgCovar());
}


//All postAnalysis functions

void Analyze::postAnalysis(){
	//Do nothing
}

void AnalyzeAverage::postAnalysis(){
  unsigned int i;
  Atom *atm;

  std::cout << std::fixed;
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

void AnalyzeEDA::postAnalysis(){
  //Take average
  if (getNData() > 1){
    getAvgCovar()/=(getNData()-1);
  }
  if (getOutput().length() > 0){
    std::ofstream out;
    out.open(getOutput().c_str(), std::ios::out);
    if (out.is_open()){
      std::cerr << "Un-diagonalized covariance matrix written to \"" << getOutput() << "\"" << std::endl;
      out << getAvgCovar() << std::endl;
      out.close();
    }
  }
  //Diagonalize
  this->diagonalizeCovar();
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

	unsigned int i,j;
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
  unsigned int i,j;
  Atom *atm;
  std::vector<Vector> coor;

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

void Analyze::averageCovariance (Molecule* cmpmol, Molecule* refmol, Eigen::MatrixXd& covarin, int& ndataIO){
  unsigned int i,j,k;
  Atom *atm;
  Eigen::MatrixXd tCovar; //Covariance at some time, t

  tCovar=Eigen::MatrixXd::Zero(covarin.rows(), covarin.cols());

  //Check selection sizes 
  if (cmpmol->getNAtomSelected() != refmol->getNAtomSelected()){
    std::cerr << std::endl << "Error: Atom number mismatch in RMSD calculation" << std::endl;
    std::cerr << "CMP-NATOM: " << cmpmol->getNAtomSelected() << ", ";
    std::cerr << "REF-NATOM: " << refmol->getNAtomSelected() << std::endl;
  }
  else{
    j=0;
    for (i=0; i< cmpmol->getNAtom(); i++){
      atm=cmpmol->getAtom(i);
      if (atm->getSel() == false){
        continue;
      }
      //Store difference (x-xavg) along diagonal
      tCovar(j,j)=atm->getX();
      j++;
      tCovar(j,j)=atm->getY();
      j++;
      tCovar(j,j)=atm->getZ();
      j++;
    }

    k=0;
    for (i=0; i< refmol->getNAtom(); i++){
      atm=refmol->getAtom(i);
      if (atm->getSel() == false){
        continue;
      }
      //Store difference (x-xavg) along diagonal
      tCovar(k,k)=tCovar(k,k)-atm->getX();
      k++;
      tCovar(k,k)=tCovar(k,k)-atm->getY();
      k++;
      tCovar(k,k)=tCovar(k,k)-atm->getZ();
      k++;
    }

    if (j != k && j != covarin.rows()){
      //No harm done, frame/structure is skipped
      std::cerr << "Warning: Atom Number Mismatch in covariance matrix!" << std::endl;
      return;
    }

    //Generate covariance matrix at time, t, for all off-diagonal elements
    //and zero diagonal afterwards
    for (i=0; i< tCovar.rows(); i++){
      for (j=i+1; j< tCovar.cols(); j++){
        tCovar(i,j)=tCovar(i,i)*tCovar(j,j);
        tCovar(j,i)=tCovar(i,j);
      }
      tCovar(i,i)=0.0;
    }

    //Add to average covariance matrix covarin="avgCovar"
    covarin+=tCovar;
    //Check possible overflow, maybe find matrix max?
      
    ndataIO++;
  }
  
}

void Analyze::projectModes(Molecule* cmpmol, Molecule* refmol, Eigen::MatrixXd& covarin){

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

void Analyze::allAnglesDihedrals(Molecule *mol, std::map<Atom*, std::vector<double> >& angin){
	unsigned int i, j;
	Chain *c;
	Atom *atmI, *iMinusTwo, *iMinusOne, *iPlusOne, *iPlusTwo, *iPlusThree;
	unsigned int size;
	std::vector<double> angles;

	for (i=0; i< mol->getChnVecSize(); i++){
    c=mol->getChain(i);
    size=c->getAtmVecSize();
    for (j=0; j< c->getAtmVecSize(); j++){
      atmI=c->getAtom(j);
			iMinusTwo=NULL;
			iMinusOne=NULL;
      iPlusOne=NULL;
      iPlusTwo=NULL;
      iPlusThree=NULL;

			if (j > 1 && atmI->getResId()-2 == c->getAtom(j-2)->getResId()){
        iMinusTwo=c->getAtom(j-2);
      }
      else{
        if (j > 1 && atmI->getResId() == c->getAtom(j-2)->getResId() && (atmI->getICode().compare(0,1,c->getAtom(j-1)->getICode(),0,1) != 0)){
          iMinusTwo=c->getAtom(j-2);
        }
      }

			if (j > 0 && atmI->getResId()-1 == c->getAtom(j-1)->getResId()){
        iMinusOne=c->getAtom(j-1);
      }
      else{
        if (j > 0 && atmI->getResId() == c->getAtom(j-1)->getResId() && (atmI->getICode().compare(0,1,c->getAtom(j-1)->getICode(),0,1) != 0)){
          iMinusOne=c->getAtom(j-1);
        }
      }

			if (j+1 < size && atmI->getResId()+1 == c->getAtom(j+1)->getResId()){
				iPlusOne=c->getAtom(j+1);
			}
			else{
				if (j+1 < size && atmI->getResId() == c->getAtom(j+1)->getResId() && (atmI->getICode().compare(0,1,c->getAtom(j+1)->getICode(),0,1) != 0)){
					iPlusOne=c->getAtom(j+1);
				}
			}

			if (j+2 < size && atmI->getResId()+2 == c->getAtom(j+2)->getResId()){
				iPlusTwo=c->getAtom(j+2);
			}
		 	else{
        if (j+2 < size && atmI->getResId() == c->getAtom(j+2)->getResId() && (atmI->getICode().compare(0,1,c->getAtom(j+2)->getICode(),0,1) != 0)){
          iPlusTwo=c->getAtom(j+2);
        }
      }

			if (j+3 < size && atmI->getResId()+3 == c->getAtom(j+3)->getResId()){
        iPlusThree=c->getAtom(j+3);
      }
			else{
        if (j+3 < size && atmI->getResId() == c->getAtom(j+3)->getResId() && (atmI->getICode().compare(0,1,c->getAtom(j+3)->getICode(),0,1) !=0)){
          iPlusThree=c->getAtom(j+3);
        }
      }
			
			angles.clear();
			angles.resize(3, 9999.9);
			if (iMinusOne != NULL && iPlusOne != NULL){
				//Get Angle
				angles.at(0)=Analyze::angle(iMinusOne->getCoor(), atmI->getCoor(), iPlusOne->getCoor());
			}
			if (iPlusOne != NULL && iPlusTwo != NULL && iPlusThree != NULL){
				//Get Dihedral
				angles.at(1)=Analyze::dihedral(atmI->getCoor(), iPlusOne->getCoor(), iPlusTwo->getCoor(), iPlusThree->getCoor());
			}
			if (iMinusTwo != NULL && iPlusTwo != NULL){
				//Get Wide Angle (i-2, i, i+2)
				angles.at(2)=Analyze::angle(iMinusTwo->getCoor(), atmI->getCoor(), iPlusTwo->getCoor());
			}
			angin.insert(std::make_pair(atmI, angles));
    }
  }
}

void Analyze::pcasso(Molecule* mol){
	Chain *c;
	Atom *ai, *aj;
	unsigned int j, start;
	double defVal;
	unsigned int icapc; //Counter for Ca or pseudocenter (pc)
	std::vector<double> iMinus6;//For non-local contacts
	std::vector<double> iPlus6; //For non-local contacts
	int diffResId;
	std::map<std::pair<Atom*, Atom*>, double> caPairDist; //Ca-Ca Distances
	std::map<std::pair<Atom*, Atom*>, double> pcPairDist; //Pseudocenter-Pseudocenter Distances
	std::map<Atom*, std::vector<double> > caAngles; //Ca-Ca Angle/Diehdral
	std::map<Atom*, std::vector<double> > pcAngles; //Pseudocenter-Pseudocenter Angle/Dihedral
	Molecule* cmol;
	std::vector<double> last;
	std::vector<double> curr;

	defVal=9999.9;
	cmol=NULL;
	last.clear();
	curr.clear();

	mol->storeSel();
	mol->select(":.CA");
	cmol=mol->clone(true,true); //Copy selection, keep original
	mol->recallSel(); //Restore original selection
	mol->eraseSel();

	//Analyze all C-alpha first
	Analyze::pairwiseDistance(cmol, caPairDist);
	Analyze::allAnglesDihedrals(cmol, caAngles);

	//Analyze all pseudocenter
	cmol->modPseudoCenter();
	Analyze::pairwiseDistance(cmol, pcPairDist);
	Analyze::allAnglesDihedrals(cmol, pcAngles);

	for (unsigned int ichain=0; ichain < cmol->getChnVecSize(); ichain++){
		c=cmol->getChain(ichain);
		for (unsigned int iatom=0; iatom < c->getAtmVecSize(); iatom++){
			ai=c->getAtom(iatom);
			//i-2, i-1, i+1, i+2, i+3, i+4, i+5 Distances
			for (icapc=0; icapc<=1; icapc++){
				//Deal with unsigned int subtraction from zero
				if (iatom == 0){
					//std::cout << defVal << " " << defVal << " ";
					curr.push_back(defVal);
					curr.push_back(defVal);
					start=iatom-0;
				}
				else if (iatom == 1){
					//std::cout << defVal << " ";
					curr.push_back(defVal);
					start=iatom-1;
				}
				else{
					start=iatom-2;
				}
				for (j=start; j<= iatom+5; j++){
					aj=c->getAtom(j);
					if (aj == NULL){
						//std::cout << defVal << " ";
						curr.push_back(defVal);
					}
					else if (j == iatom){
						//Distance == 0
						continue;
					}
					else{
						if (icapc == 0){ //Ca
							//std::cout << caPairDist.at(std::make_pair(ai, aj)) << " ";
							curr.push_back(caPairDist.at(std::make_pair(ai, aj)));
						}
						else{ //Pseudocenter
							//std::cout << pcPairDist.at(std::make_pair(ai, aj)) << " ";
							curr.push_back(pcPairDist.at(std::make_pair(ai, aj)));
						}
					}
				}
			}

			//Angles and Dihedrals
			for (icapc=0; icapc<=1; icapc++){
				if (icapc == 0){ //Ca
					for (j=0; j< caAngles.at(ai).size(); j++){
						//std::cout << caAngles.at(ai).at(j) << " ";
						curr.push_back(caAngles.at(ai).at(j));
					}
				}
				else{ //Pseudocenter
					for (j=0; j< pcAngles.at(ai).size(); j++){
            //std::cout << pcAngles.at(ai).at(j) << " ";
						curr.push_back(pcAngles.at(ai).at(j));
          }
				}
			}
		
			//Shortest non-local contact distance, >= i+6 and <= i-6
			iPlus6.clear();
			iMinus6.clear();
			for (j=0; j< cmol->getAtmVecSize(); j++){
				aj=cmol->getAtom(j);
				if (ai == aj){
					continue;
				}
				//i+6
				if (ai->getChainId().compare(aj->getChainId()) != 0){
					//atom i and atom j are on different chains
					iPlus6.push_back(caPairDist.at(std::make_pair(ai, aj)));
				}
				else{
					//atom i and atom j are on the same chain
					diffResId=aj->getResId() - ai->getResId();
					if (diffResId >= 6){
						iPlus6.push_back(caPairDist.at(std::make_pair(ai, aj)));
					}
					else{
						if (diffResId == 0 && (Misc::atoi(aj->getICode()) - Misc::atoi(ai->getICode()) >= 6)){
							iPlus6.push_back(caPairDist.at(std::make_pair(ai, aj)));
						}
					}
				}

				//i-6
				if (ai->getChainId().compare(aj->getChainId()) != 0){
					//atom i and atom j are on different chains
					iMinus6.push_back(caPairDist.at(std::make_pair(ai, aj)));
				}
				else{
					//atom i and atom j are on the same chain
					diffResId=aj->getResId() - ai->getResId();
					if (diffResId <= -6){
						iMinus6.push_back(caPairDist.at(std::make_pair(ai, aj)));
					}
					else{
						if (diffResId == 0 && (Misc::atoi(aj->getICode()) - Misc::atoi(ai->getICode()) <= -6)){
							iMinus6.push_back(caPairDist.at(std::make_pair(ai, aj)));
						}
					}
				}
			}
			std::sort(iPlus6.begin(),iPlus6.end());
			for (j=0; j< 3; j++){
				if (j < iPlus6.size()){
					//std::cout << iPlus6.at(j) << " ";
					curr.push_back(iPlus6.at(j));
				}
				else{
					//std::cout << defVal << " ";
					curr.push_back(defVal);
				}
			}
			std::sort(iMinus6.begin(),iMinus6.end());
			for (j=0; j< 3; j++){
        if (j < iMinus6.size()){
          //std::cout << iMinus6.at(j) << " ";
					curr.push_back(iMinus6.at(j));
        }
        else{
          //std::cout << defVal << " ";
					curr.push_back(defVal);
        }
      }

			//Print output
			if (iatom > 0){
				//Not first atom of chain
				//Print S(i+1) for last atom
				for (j=0; j< curr.size(); j++){
          std::cout << curr.at(j) << " ";
        }
				std::cout << std::endl;
				
				std::cout << ai->getPdbId() << " " << ai->getChainId() << " " << ai->getResId() << " ";
				//Print S(i)
				for (j=0; j< curr.size(); j++){
          std::cout << curr.at(j) << " ";
        }
				//Print S(i-1) stored in std::vector "last"
				for (j=0; j< last.size(); j++){
					std::cout << last.at(j) << " ";
				}
				if (iatom+1 == c->getAtmVecSize()){
					//Last atom of chain
					//Print S(i+1) which is all default values
					for (j=0; j< curr.size(); j++){
	          std::cout << defVal << " ";
	        }
					std::cout << std::endl;
				}
			}
			else{
				//First atom in chain, last is empty
				std::cout << ai->getPdbId() << " " << ai->getChainId() << " " << ai->getResId() << " ";
				//Print S(i)
				for (j=0; j< curr.size(); j++){
          std::cout << curr.at(j) << " ";
        }
				//Print S(i-1) which is all default values
				for (j=0; j< curr.size(); j++){
          std::cout << defVal << " ";
        }
			}
			last=curr;
			curr.clear();

		}
	}

	if (cmol != NULL){
		delete cmol;
	}
}
