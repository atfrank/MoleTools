//Sean M. Law
//Aaron T. Frank
    
/*
This file is part of MoleTools.

MoleTools is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MoleTools is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MoleTools.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "Analyze.hpp"

#include "Molecule.hpp"
#include "Chain.hpp"
#include "Residue.hpp"
#include "Atom.hpp"
#include "Coor.hpp"
#include "Constants.hpp"
#include "DTree.hpp"
#include "PCASSO.hpp"
#include "Misc.hpp"

#include <fstream>
#include <iomanip>

Analyze::Analyze (){
  sel.clear();
  mol.clear();
  tdata.clear();
  fdata.clear();
  ndata=0;
  resel=false;
  avgCovar.resize(0,0);
  ifile.clear();
  ofile.clear();
  verbose=false;
}

Analyze::~Analyze (){
  //Do nothing
}

AnalyzePcasso::AnalyzePcasso (std::string delim){
  std::vector<std::string> tokens;

  pout=PREDICT;
  t.clear();
  t.resize(PCASSO::getNTree());
  for (unsigned int i=0; i< PCASSO::getNTree(); i++){
    t.at(i)=new DTree;
    Misc::splitStr(Misc::trim(PCASSO::getTree(i)), " \t", tokens, false);
    t.at(i)->genDTree(tokens, delim);
  }
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

void Analyze::setMol(const int& element, Molecule* molin){
  this->mol.at(element)=molin;
}

void Analyze::clearMol(){
  this->mol.clear();
}

void Analyze::resizeNMol(const int sizein){
  this->mol.resize(sizein);
}

void Analyze::readTopology(Molecule* molin, std::string topin){
  if (topin.length() > 0){
    molin->readTopology(topin);
    molin->setMass();
    molin->setCharge();
  }
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

std::vector<std::vector<double> >& Analyze::getFDataVec(){
  return fdata;
}

Eigen::MatrixXd& Analyze::getCovar(){
  return avgCovar;
}

void Analyze::addModes(const std::vector<unsigned int>& modesin){
  modes=modesin;
}

std::vector<unsigned int>& Analyze::getModes(){
  return modes;
}

void Analyze::initCovar(const unsigned int &xin, const unsigned int &yin){
  avgCovar=Eigen::MatrixXd::Zero(xin, yin);
}

Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> Analyze::diagonalizeCovar(){
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
  
  return eigOut;
}

void Analyze::setEigenMode(const unsigned int& modein){
  //Sets the molecule to the first mode
  double eigval;
  Eigen::MatrixXd eigvec;
  unsigned int nrow;
  unsigned int i, j;
  Atom* a;

  if (this->getNMol() > 1){
    this->setMol(1, this->getMol(0)->clone());
  }
  else{
    this->setupMolSel(this->getMol(0)->clone());
  }
  
  nrow=this->getEigen().eigenvalues().rows();
  eigval=this->getEigen().eigenvalues()[nrow-modein];
  //Note that the amplitude is SQRT(Eigenvalue) since the eigenvalue
  //is a measure of the variance while the amplitude is a measure of
  //the standard deviation!
  eigvec=this->getEigen().eigenvectors().col(nrow-modein)*sqrt(eigval);
  if (nrow != 3*this->getMol(1)->getAtmVecSize()){
    std::cerr << "Warning: The number of 3N coordinates (" << 3*this->getMol(1)->getAtmVecSize();
    std::cerr << ") and the number of elements in the covariance matrix (" << nrow << ") do not match!" << std::endl;
  }
  else{
    i=0;
    for (j=0; j< this->getMol(1)->getAtmVecSize(); j++){
      a=this->getMol(1)->getAtom(j);
      a->setX(a->getX()+eigvec(i));
      i++;
      a->setY(a->getY()+eigvec(i));
      i++;
      a->setZ(a->getZ()+eigvec(i));
      i++;
    }
  }
}

void Analyze::writeEigenOverlap(Analyze* cmpin, std::vector<unsigned int>& modein){
  unsigned int i, j;
  Eigen::MatrixXd refvec;
  Eigen::MatrixXd cmpvec;
  double overlap;

  refvec=this->getEigen().eigenvectors();
  cmpvec=cmpin->getEigen().eigenvectors();


  if (modein.size() == 0 || (modein.size() == 1 && modein.at(0) == 0)){
    //If mode = 0 then print all eigenvector overlaps
    modein.clear();
    for (j=0; j< static_cast<unsigned int>(refvec.cols()); j++){
       modein.push_back(j+1);
    }
  }

  std::cout << "#" << this->getInput() << "  " << cmpin->getInput() << "  " << "-1 <= Overlap <= 1" << std::endl;
  if (refvec.cols() > 0 && refvec.cols() == cmpvec.cols()){
    if (refvec.col(1).rows() == cmpvec.col(1).rows()){
      for (i=0; i< static_cast<unsigned int>(refvec.cols()); i++){
        if (std::find(modein.begin(), modein.end(), i+1) == modein.end()){
          continue;
        }
        for (j=0; j< static_cast<unsigned int>(cmpvec.cols()); j++){
          if (std::find(modein.begin(), modein.end(), j+1) == modein.end()){
            continue;
          }
          overlap=refvec.col(refvec.cols()-(i+1)).dot(cmpvec.col(cmpvec.cols()-(j+1)));
          std::cout << i+1 << "  " << j+1 << "  " << overlap << std::endl;
        }
        std::cout << std::endl;
      }
    }
  }
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

void Analyze::setVerbose(bool verbosein){
  verbose=verbosein;
}

bool Analyze::getVerbose(){
  return verbose;
}

void Analyze::setEigen(const Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd>& eigenin){
  eigen=eigenin;
}

Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd>& Analyze::getEigen(){
  return eigen;
}

void Analyze::readCovariance(){
  std::ifstream inpFile;
  std::istream* inp;
  std::string line;
  unsigned int nrow;
  std::vector<double> tmp;

  //Read covariance matrix
  if (getInput().length() > 0){
    nrow=0;

    inpFile.open(getInput().c_str(), std::ios::in);
    inp=&inpFile;

    while (inp->good() && !(inp->eof())){
      getline(*inp, line);
      if (line.length() > 0){
        Misc::splitNum(Misc::trim(line), " \t", tmp, false);
        if (nrow == 0){
          this->initCovar(tmp.size(), tmp.size());
        }
        for (unsigned int ncol=0; ncol< tmp.size(); ncol++){
          getCovar()(nrow, ncol)=tmp.at(ncol);
        }
        nrow++;
      }
    }
    if (inpFile.is_open()){
      inpFile.close();
    }
    if (nrow == 0){
      std::cerr << std::endl << "Warning: The covariance matrix could not be read" << std::endl;
    }
  }
  else{
    std::cerr << "Warning: A covariance matrix was not provided" << std::endl;
  }
}




//All preAnalysis Functions
void Analyze::preAnalysis(){
  //Do nothing
}

void Analyze::preAnalysis(Molecule* molin, std::string topin){
  this->readTopology(molin, topin);
  this->setupMolSel(molin);
}

void AnalyzeRMSD::preAnalysis(Molecule* molin, std::string topin){
  this->readTopology(molin, topin);
  this->setupMolSel(molin);
  Molecule* refmol=molin->clone();
  this->readTopology(refmol, topin);
  this->setupMolSel(refmol);
}

void AnalyzeRMSF::preAnalysis(Molecule* molin, std::string topin){
  this->readTopology(molin, topin);
  this->setupMolSel(molin);
  Molecule* refmol=molin->clone();
  this->readTopology(refmol, topin);
  this->setupMolSel(refmol);
}

void AnalyzeAverage::preAnalysis(Molecule* molin, std::string topin){
  this->readTopology(molin, topin);
  this->setupMolSel(molin);
  Molecule* refmol=molin->clone();
  this->readTopology(refmol, topin);
  this->setupMolSel(refmol);
  //Zero coordinates for future appending
  this->getMol(1)->zeroCoor();
}

void AnalyzeCovariance::preAnalysis(Molecule* molin, std::string topin){
  this->readTopology(molin, topin);
  this->setupMolSel(molin);
  Molecule* refmol=molin->clone();
  this->readTopology(refmol, topin);
  this->setupMolSel(refmol);

  if (getInput().length() != 0){
    readCovariance();

    //Diagonalize covariance matrix
    setEigen(diagonalizeCovar());
  }
  else{
    //Resize matrix to 3N x 3N and zero
    this->initCovar(3*refmol->getNAtomSelected(), 3*refmol->getNAtomSelected());
  }
}

void AnalyzeCovariance::preAnalysis(){

  readCovariance();

  //Diagonalize covariance matrix
  setEigen(diagonalizeCovar());
}

void AnalyzeProjection::preAnalysis(Molecule* molin, std::string topin){
  std::ifstream inpFile;
  std::istream* inp;
  std::string line;
  unsigned int nrow;
  std::vector<double> tmp;
  unsigned int N3;

  this->readTopology(molin, topin);
  this->setupMolSel(molin);
  Molecule* refmol=molin->clone();
  this->readTopology(refmol, topin);
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
          getCovar()(nrow, ncol)=tmp.at(ncol);
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

  //Diagonalize covariance matrix
  setEigen(diagonalizeCovar());

  //std::cout << getEigen().eigenvalues() << std::endl;
  //std::cout << getEigen().eigenvectors().col(N3-1) << std::endl;
}

void AnalyzePcasso::preAnalysis(Molecule* molin, std::string fin){
  this->setupMolSel(molin);
  this->setInput(fin);
  
  //Resize FData in Analyze::pcasso
}


//All runAnalysis Functions

void AnalyzeCOG::runAnalysis(){
  Coor xyz=Analyze::centerOfGeometry(this->getMol(0));
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

void AnalyzeCovariance::runAnalysis(){
  Analyze::averageCovariance(this->getMol(0), this->getMol(1), getCovar(), getNData());
}

void AnalyzeProjection::runAnalysis(){
  std::vector<double> pc;

  pc=Analyze::projectModes(this->getMol(0), this->getMol(1), getEigen(), getModes());

  for (unsigned int i=0; i< pc.size(); i++){
    std::cout << std::fixed;
    std::cout << std::setw(9) << std::right << std::setprecision(3) << pc.at(i);
  }
}

void AnalyzeGyrationTensor::runAnalysis(){
  Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eigVal;
  unsigned int nrow;

  eigVal.compute(Analyze::gyrationTensor(this->getMol(0)));
  nrow=eigVal.eigenvalues().rows();

  for (unsigned int i=1; i<= 3; i++){
    std::cout << std::fixed;
    std::cout << std::setw(9) << std::right << std::setprecision(3) << sqrt(eigVal.eigenvalues()[nrow-i]);
  }
}

void AnalyzeRadiusOfGyration::runAnalysis(){
  Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eigVal;

  eigVal.compute(Analyze::gyrationTensor(this->getMol(0)));

  std::cout << std::fixed;
  std::cout << std::setw(9) << std::right << std::setprecision(3) << sqrt(eigVal.eigenvalues()[0]+eigVal.eigenvalues()[1]+eigVal.eigenvalues()[2]);
  //Rgyr = sqrt(Lambda1 + Lambda2 +Lambda3)

  std::cerr << std::endl << "Warning: The radius of gyration is not mass-weighted!" << std::endl;
}

void AnalyzeEllipsoid::runAnalysis(){
  Coor xyz;
  Molecule* cogmol;
  Chain* c;
  Residue* r;
  Atom* a;
  Eigen::Matrix3d tensor;
   
  tensor=Analyze::gyrationTensor(this->getMol(0));
 
  cogmol=new Molecule;
  c=new Chain;
  r=new Residue;
  a=new Atom;

  xyz=Analyze::centerOfGeometry(this->getMol(0));

  a->dummy();
  a->setCoor(xyz);
  r->addAtom(a);
  c->addAtom(a);
  c->addResidue(r);
  cogmol->addAtom(a);
  cogmol->addResidue(r);
  cogmol->addChain(c);

  //Outputs in the order of U[0,0], U[1,1], U[2,2], U[0,1], U[0,2], U[1,2]
  std::cout << std::fixed;
  std::cout << std::setw(8) << std::right << std::setprecision(0) << tensor(0,0)*1E4;
  std::cout << std::setw(8) << std::right << std::setprecision(0) << tensor(1,1)*1E4;
  std::cout << std::setw(8) << std::right << std::setprecision(0) << tensor(2,2)*1E4;
  std::cout << std::setw(8) << std::right << std::setprecision(0) << tensor(0,1)*1E4;
  std::cout << std::setw(8) << std::right << std::setprecision(0) << tensor(0,2)*1E4;
  std::cout << std::setw(8) << std::right << std::setprecision(0) << tensor(1,2)*1E4;

  delete cogmol;
}

void AnalyzePairwiseDistance::runAnalysis(){
//  std::map<std::pair<Atom*, Atom*>, double> pdist;
  std::vector<std::vector<double> > pdist;
  Atom* ai;
  Atom* aj;
  
  this->getMol(0)->assignAtmInx();
  Analyze::pairwiseDistance(this->getMol(0), pdist);

  for (unsigned int i=0; i< this->getMol(0)->getAtmVecSize(); i++){
    ai=this->getMol(0)->getAtom(i);
    for (unsigned int j=i+1; j< this->getMol(0)->getAtmVecSize(); j++){
      aj=this->getMol(0)->getAtom(j);
//      std::cout << "  " << ai->getSummary() << "-" << aj->getSummary();
      std::cout << std::fixed;
      std::cout << std::setw(9) << std::right << std::setprecision(3) << pdist.at(ai->getAtmInx()).at(aj->getAtmInx());
    }
  }
}

void AnalyzePcasso::runAnalysis(){
  std::ifstream dsspFile;
  std::istream* dsspinp;
  std::string line;
  std::vector<std::string> dssp;
  std::vector<std::string> s;
  unsigned int natom;
  bool last;

  Analyze::pcasso(this->getMol(0), this->getFDataVec()); //PCASSO features get stored in the second argument (a 2-D vector double)

  natom=0;
  last=false;

  //Read DSSP file first
  if (this->getInput().length() > 0){
    dsspFile.open(this->getInput().c_str(), std::ios::in);
    dsspinp=&dsspFile;
    while (dsspinp->good()){
      getline(*dsspinp, line);
      if (line.length() == 0){
        continue;
      }
      Misc::splitStr(line, " \t", s, false);
      if (s.size() == 2 && last==false){
        if (s.at(s.size()-2).compare(0,1,"E") == 0 || s.at(s.size()-2).compare(0,1,"B") == 0){
          dssp.push_back("E");
        }
        else if (s.at(s.size()-2).compare(0,1,"H") == 0 || s.at(s.size()-2).compare(0,1,"I") == 0 || s.at(s.size()-2).compare(0,1,"G") == 0){
          dssp.push_back("H");
        }
        else if (s.at(s.size()-2).compare(0,1,"-") == 0 || s.at(s.size()-2).compare(0,1,"S") == 0 || s.at(s.size()-2).compare(0,1,"T") == 0){
          dssp.push_back("C");
        }
        else{
          std::cerr << "Warning unrecognized DSSP classification \"" << s.at(s.size()-2) << "\" has been set to \"X\"" << std::endl;
          dssp.push_back("X");
        }
      }
      else{
        last=true;
        if (s.size() > 0){
          dssp.push_back(s.at(s.size()-1));
        }
        else{
          dssp.push_back("");
        } 
      }
    }
    if (dsspFile.is_open()){
      dsspFile.close();
    }
  }

  std::vector<std::vector<double> > &feat=this->getFDataVec();
  std::map<std::string, unsigned int> vote; 
  std::map<std::string, unsigned int>::iterator iter;
  std::string tmpClass;
  std::string maxClass;
  unsigned int maxVote;
  bool majority;
  unsigned int ntree;

  ntree=PCASSO::getNTree();

  if (this->getOutType() == PREDICT || this->getOutType() == PREDICTION){
    for (unsigned int i=0; i< feat.size(); i++){
      vote.clear();
      maxVote=0;
      majority=false;
      for (unsigned int j=0; j< ntree && majority == false; j++){
        tmpClass=t.at(j)->getDTreeClass(feat.at(i));
        if (vote.find(tmpClass) != vote.end()){
          vote.at(tmpClass)++;
          if (vote.at(tmpClass) > maxVote){
            //Find majority vote, method adapted from openCV
            //Is pseudo-random since it depends on the order of the trees
            maxVote=vote.at(tmpClass);
            maxClass=tmpClass;
            if (static_cast<float>(maxVote)/ntree > 0.5){ //Unsigned integer division!
              majority=true;
            }
          }
        }
        else{
          vote.insert(std::pair<std::string, unsigned int>(tmpClass,1));
        }
      }
      if (this->getVerbose() == true && this->getMol(0)->getNAtom() == feat.size()){
        std::cout << " " << this->getMol(0)->getAtom(i)->getSummary();
        std::cout << " " << maxClass ; //Print majority vote
      }
      else{
        std::cout << " " << maxClass; //Print majority vote
      }
    }
  }
  else if (this->getOutType() == FEATURES || this->getOutType() == FEATURE){
    //Print features
    for (unsigned int i=0; i < feat.size(); i++){
      if (this->getInput().length() > 0 && natom < dssp.size()){
        std::cout << dssp.at(natom) << ",";
      }
      else{
        std::cout << "X" << ",";
      }

      for (unsigned int j=0; j< feat.at(i).size(); j++){
        std::cout << feat.at(i).at(j);
        if (j < feat.at(i).size()-1){
          std::cout << ",";
        }
      }

      std::cout << std::endl;
      natom++;
    } //Loop through chains

    if (this->getInput().length() > 0 && dssp.size() != natom){
      std::cerr << "Warning: DSSP (" << dssp.size() << ") and NATOM (" << natom << ") mismatch" << std::endl;
    }
  }
  else{
    std::cerr << "Warning: Unrecognized PCASSO output type" << std::endl;
  }

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
    tdata.at(j)=sqrt(tdata.at(j)/(getNData()-1));
    std::cout << atm->getSummary();
    std::cout << std::setw(9) << std::right << std::setprecision(3) << tdata.at(j);
    std::cout << std::endl;
    j++;
  }
}

void AnalyzeCovariance::postAnalysis(){
  //Take average
  if (getNData() > 1){
    getCovar()/=(getNData()-1);
  }
  if (getOutput().length() > 0){
    std::ofstream out;
    out.open(getOutput().c_str(), std::ios::out);
    if (out.is_open()){
      std::cerr << "Un-diagonalized covariance matrix written to \"" << getOutput() << "\"" << std::endl;
      out << getCovar() << std::endl;
      out.close();
    }
  }
}

void AnalyzePcasso::postAnalysis(){
  while (!t.empty()){
    t.back()->delDTree();
    t.pop_back();
  }
}

//Basic analysis functions

Coor Analyze::centerOfGeometry(Molecule *mol, bool selFlag){
  Coor cog=Coor(0.0, 0.0, 0.0);

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
  std::vector<Coor> coor;
  Coor d;

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
  std::vector<Coor> coor;

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

double Analyze::distance (const Coor& u, const Coor& v){
  Coor d=u-v;
  return d.norm();
}

double Analyze::angle (const Coor& u, const Coor& v, const Coor& w){
  double angle;
  Coor dx, dy;
  double dp, nx, ny;

  dx=u-v;
  dy=w-v;

  nx=dx.norm();
  ny=dy.norm();

  dp=dx.dot(dy);

  angle=acos(dp/(nx*ny));

  return angle/PI*180.0;
}

double Analyze::dihedral (const Coor& t, const Coor& u, const Coor& v, const Coor& w) {
  double dihedral;
  Coor dx, dy, dz, p1, p2, p3;
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

    if (j != k && j != static_cast<unsigned int>(covarin.rows())){
      //No harm done, frame/structure is skipped
      std::cerr << "Warning: Atom Number Mismatch in covariance matrix!" << std::endl;
      return;
    }

    //Generate covariance matrix at time, t, for all off-diagonal elements
    //Note that the diganonal is NOT zero!
    for (i=0; i< static_cast<unsigned int>(tCovar.rows()); i++){
      for (j=i+1; j< static_cast<unsigned int>(tCovar.cols()); j++){
        tCovar(i,j)=tCovar(i,i)*tCovar(j,j);
        tCovar(j,i)=tCovar(i,j);
      }
      tCovar(i,i)=tCovar(i,i)*tCovar(i,i);
    }

    //Add to average covariance matrix covarin="avgCovar"
    covarin+=tCovar;
    //Check possible overflow, maybe find matrix max?
      
    ndataIO++;
  }
  
}

std::vector<double> Analyze::projectModes(Molecule* cmpmol, Molecule* refmol, const Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd>& eigenin, const std::vector<unsigned int>& modesin){
  unsigned int i,j,k,m;
  Atom *atm;
  std::vector<double> proj;
  double defVal;
  Eigen::VectorXd dt; //Difference, (x(t)-xavg), at some time, t

  dt=Eigen::VectorXd::Zero(eigenin.eigenvalues().rows());
  proj.clear();
  proj.resize(modesin.size(),0);
  defVal=9999.9;

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
      dt(j)=atm->getX();
      j++;
      dt(j)=atm->getY();
      j++;
      dt(j)=atm->getZ();
      j++;
    }

    k=0;
    for (i=0; i< refmol->getNAtom(); i++){
      atm=refmol->getAtom(i);
      if (atm->getSel() == false){
        continue;
      }
      //Store difference (x-xavg) along diagonal
      dt(k)=dt(k)-atm->getX();
      k++;
      dt(k)=dt(k)-atm->getY();
      k++;
      dt(k)=dt(k)-atm->getZ();
      k++;
    }

    if (j != k && j != static_cast<unsigned int>(eigenin.eigenvalues().rows())){
      //No harm done, frame/structure is skipped
      std::cerr << "Warning: Atom Number Mismatch in covariance matrix!" << std::endl;
      return proj;
    }
  
    //Perform dot product for each requested mode
    for (m=0; m< modesin.size(); m++){
      if (m >= static_cast<unsigned int>(eigenin.eigenvalues().rows())){
        std::cerr << "Warning: Mode " << modesin.at(m) << " does not exist and projection was set to default " << defVal << std::endl;
        proj.at(m)=defVal;
      }
      else{
        //Note that mode 1 is the last eigenvector and modesin uses base one
        proj.at(m)=eigenin.eigenvectors().col(eigenin.eigenvalues().rows()-modesin.at(m)).dot(dt);
      }
    }
  }
  
  return proj;
}

void Analyze::pairwiseDistance(Molecule *mol, std::vector<std::vector<double> >& pdin){
  std::vector<Atom*>::iterator ai;
  std::vector<Atom*>::iterator aj;
  unsigned int natom;
  unsigned int aiInx;
  bool flag;

  natom=mol->getAtmVecSize();

  pdin.clear();
  pdin.resize(natom);

  for (ai=mol->getAtmVec().begin(); ai != mol->getAtmVec().end(); ++ai){
    aiInx=(*ai)->getAtmInx();
    pdin.at(aiInx).resize(natom);
    pdin.at(aiInx).at(aiInx)=0.0; //Zero diagonal
  }

  for (ai=mol->getAtmVec().begin(); ai != mol->getAtmVec().end(); ++ai){
    aiInx=(*ai)->getAtmInx();
    if ((*ai)->getX() < 9999.9){
      flag=true;
    }
    else{
      flag=false;
    }
    //Lower Triangle
    for (aj=mol->getAtmVec().begin(); aj != ai; ++aj){
      if (flag && (*aj)->getX() < 9999.9){
        pdin.at(aiInx).at((*aj)->getAtmInx())=Analyze::distance((*ai)->getCoor(), (*aj)->getCoor());
      }
      else{
        pdin.at(aiInx).at((*aj)->getAtmInx())=9999.9;
      }
    }
  
    //Upper Triangle
    for (aj=ai+1; aj != mol->getAtmVec().end(); ++aj){
      if (flag && (*aj)->getX() < 9999.9){
        pdin.at(aiInx).at((*aj)->getAtmInx())=Analyze::distance((*ai)->getCoor(), (*aj)->getCoor());
      }
      else{
        pdin.at(aiInx).at((*aj)->getAtmInx())=9999.9;
      }
    }
  }

}

void Analyze::allAnglesDihedrals(Molecule *mol, std::vector<std::vector<double> >& anglesin){
  unsigned int i, j;
  Chain *c;
  Atom *atmI, *iMinusTwo, *iMinusOne, *iPlusOne, *iPlusTwo, *iPlusThree;
  unsigned int size;
  std::vector<double> angles;

  anglesin.resize(mol->getAtmVecSize());

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
      anglesin.at(atmI->getAtmInx()).resize(3);
      anglesin.at(atmI->getAtmInx())=angles;
    }
  }
}

void Analyze::pcasso(Molecule* mol, std::vector<std::vector<double> > &fdataIO){
  Chain *c;
  Atom *ai, *aj, *ak;
  unsigned int i, j, start;
  double defVal;
  double iMinus6;//For non-local contacts
  double iPlus6; //For non-local contacts
  unsigned minx;
  unsigned pinx;
  int diffResId;
  std::vector<std::vector<double> > caPairDist; //Ca-Ca Distances
  std::vector<std::vector<double> > pcPairDist; //Pc-Pc Distances
  std::vector<std::vector<double> > caAngles; //Ca-Ca Angle/Diehdral/Wide
  std::vector<std::vector<double> > pcAngles; //Pc-Pc Angle/Dihedral/Wide
  Molecule *camol, *pcmol;
  unsigned int natom;
  double dist;

  defVal=9999.9;
  camol=NULL;
  pcmol=NULL;

  mol->storeSel();
  mol->select(":.CA");
  camol=mol->clone(true,true); //Copy selection, keep original
  pcmol=mol->clone(true,true);
  mol->recallSel(); //Restore original selection
  mol->eraseSel();

  camol->assignAtmInx();
  pcmol->assignAtmInx();

  //Analyze all C-alpha first
  Analyze::pairwiseDistance(camol, caPairDist);
  Analyze::allAnglesDihedrals(camol, caAngles);

  natom=camol->getAtmVecSize();

  for (unsigned int ichain=0; ichain < camol->getChnVecSize(); ichain++){
    c=camol->getChain(ichain);
    for (unsigned int iatom=0; iatom < c->getAtmVecSize(); iatom++){
      ai=c->getAtom(iatom);
      ai->clearData();
      i=ai->getAtmInx();
      //i-5, i-4, i-3, i-2, i-1, i+1, i+2, i+3, i+4, i+5 Distances
      //Deal with unsigned int subtraction from zero
      if (iatom == 0){
        ai->addData(defVal);
        ai->addData(defVal);
        ai->addData(defVal);
        ai->addData(defVal);
        ai->addData(defVal);
        start=iatom-0;
      }
      else if (iatom == 1){
        //std::cout << defVal << " ";
        ai->addData(defVal);
        ai->addData(defVal);
        ai->addData(defVal);
        ai->addData(defVal);
        start=iatom-1;
      }
      else if (iatom == 2){
        ai->addData(defVal);
        ai->addData(defVal);
        ai->addData(defVal);
        start=iatom-2;
      }
      else if (iatom == 3){
        ai->addData(defVal);
        ai->addData(defVal);
        start=iatom-3;
      }
      else if (iatom == 4){
        ai->addData(defVal);  
        start=iatom-4;
      }
      else{
        start=iatom-5;
      }
      for (j=start; j<= iatom+5; j++){
        aj=c->getAtom(j);
        if (aj == NULL){
          ai->addData(defVal);
        }
        else if (j == iatom){
          //Distance == 0
          continue;
        }
        else{
          ai->addData(caPairDist.at(i).at(aj->getAtmInx()));
        }
      }

      //Angles and Dihedrals
      for (j=0; j< caAngles.at(i).size(); j++){
        ai->addData(caAngles.at(i).at(j));
      }
    
      //Shortest non-local contact distance, >= i+6 and <= i-6
      iPlus6=1E10;
      iMinus6=1E10;
      pinx=natom;
      minx=natom;

      for (j=0; j< natom; j++){
        aj=camol->getAtom(j);

        if (ai == aj){
          continue;
        }
        
        dist=caPairDist.at(i).at(j);

        //i+6
        //Assess distance first to avoid unnecessary string comparison
        if (dist < iPlus6){
          if (ai->getChainId().compare(aj->getChainId()) != 0){
            //atom i and atom j are on different chains
            pinx=j;
            iPlus6=dist;
          }
          else{
            //atom i and atom j are on the same chain
            diffResId=aj->getResId() - ai->getResId();
            if (diffResId >= 6){
              pinx=j;
              iPlus6=dist;
            }
            else{
              if (diffResId == 0 && (Misc::atoi(aj->getICode()) - Misc::atoi(ai->getICode()) >= 6)){
                pinx=j;
                iPlus6=dist;
              }
            }
          }
        }

        //i-6
        //Assess distance first to avoid unnecessary string comparison
        if (dist < iMinus6){
          if (ai->getChainId().compare(aj->getChainId()) != 0){
          //atom i and atom j are on different chains
            minx=j;
            iMinus6=dist;
          }
          else{
            //atom i and atom j are on the same chain
            diffResId=aj->getResId() - ai->getResId();
            if (diffResId <= -6){
              minx=j;
              iMinus6=dist;
            }
            else{
              if (diffResId == 0 && (Misc::atoi(aj->getICode()) - Misc::atoi(ai->getICode()) <= -6)){
                minx=j;
                iMinus6=dist;
              }
            }
          }
        }
      }

      int k;
      unsigned int q;
      int max;
      max=10;
      if (iatom == 0){
        for (k=0; k< max; k++){
          ai->addData(defVal);
        }
        start=0;
      }
      else{
        start=iatom-1;
      }
      for (q=start; q<= iatom+1; q++){
        if (q < c->getAtmVecSize()){
          ak=c->getAtom(q);
          for (k=static_cast<int>(pinx)-2; k<=static_cast<int>(pinx)+2; k++){
            if (k >= 0 && k < static_cast<int>(natom)){
              aj=camol->getAtom(k);
              ai->addData(caPairDist.at(ak->getAtmInx()).at(aj->getAtmInx()));
            }
            else{
              ai->addData(defVal);
            }
          }
          for (k=static_cast<int>(minx)-2; k<=static_cast<int>(minx)+2; k++){
            if (k >= 0 && k < static_cast<int>(natom)){
              aj=camol->getAtom(k);
              ai->addData(caPairDist.at(ak->getAtmInx()).at(aj->getAtmInx()));
            }
            else{
              ai->addData(defVal);
            } 
          }
        }
        else{
          for (k=0; k< max; k++){
            ai->addData(defVal);
          }
        }
      }

    } //Loop through atoms
  }//Loop through chains

  //Analyze all pseudocenter
  pcmol->modPseudoCenter();
  Analyze::pairwiseDistance(pcmol, pcPairDist);
  Analyze::allAnglesDihedrals(pcmol, pcAngles);

  natom=pcmol->getAtmVecSize();

  for (unsigned int ichain=0; ichain < pcmol->getChnVecSize(); ichain++){
    c=pcmol->getChain(ichain);
    for (unsigned int iatom=0; iatom < c->getAtmVecSize(); iatom++){
      ai=camol->getChain(ichain)->getAtom(iatom); //From C-alpha
      ak=c->getAtom(iatom);
      i=ak->getAtmInx();
      //i-5, i-4, i-3, i-2, i-1, i+1, i+2, i+3, i+4, i+5 Distances
      //Deal with unsigned int subtraction from zero
      if (iatom == 0){
        ai->addData(defVal);
        ai->addData(defVal);
        ai->addData(defVal);
        ai->addData(defVal);
        ai->addData(defVal);
        start=iatom-0;
      }
      else if (iatom == 1){
        //std::cout << defVal << " ";
        ai->addData(defVal);
        ai->addData(defVal);
        ai->addData(defVal);
        ai->addData(defVal);
        start=iatom-1;
      }
      else if (iatom == 2){
        ai->addData(defVal);
        ai->addData(defVal);
        ai->addData(defVal);
        start=iatom-2;
      }
      else if (iatom == 3){
        ai->addData(defVal);
        ai->addData(defVal);
        start=iatom-3;
      }
      else if (iatom == 4){
        ai->addData(defVal);  
        start=iatom-4;
      }
      else{
        start=iatom-5;
      }
      for (j=start; j<= iatom+5; j++){
        aj=c->getAtom(j);
        if (aj == NULL){
          ai->addData(defVal);
        }
        else if (j == iatom){
          //Distance == 0
          continue;
        }
        else{
          ai->addData(pcPairDist.at(i).at(aj->getAtmInx()));
        }
      }

      //Angles and Dihedrals
      for (j=0; j< pcAngles.at(ak->getAtmInx()).size(); j++){
        ai->addData(pcAngles.at(ak->getAtmInx()).at(j));
      }
    
      //Shortest non-local contact distance, >= i+6 and <= i-6
      iPlus6=1E10;
      iMinus6=1E10;
      pinx=natom;
      minx=natom;
      for (j=0; j< natom; j++){
        aj=pcmol->getAtom(j);
        if (ak == aj){
          continue;
        }

        dist=pcPairDist.at(i).at(j);

        //i+6
        //Assess distance first to avoid unnecessary string comparison
        if (dist < iPlus6){
          if (ak->getChainId().compare(aj->getChainId()) != 0){
            //atom i and atom j are on different chains
            pinx=j;
            iPlus6=dist;
          }
          else{
            //atom i and atom j are on the same chain
            diffResId=aj->getResId() - ak->getResId();
            if (diffResId >= 6){
              pinx=j;
              iPlus6=dist;
            }
            else{
              if (diffResId == 0 && (Misc::atoi(aj->getICode()) - Misc::atoi(ak->getICode()) >= 6)){
                pinx=j;
                iPlus6=dist;
              }
            }
          }
        }

        //i-6
        //Assess distance first to avoid unnecessary string comparison
        if (dist < iMinus6){
          if (ak->getChainId().compare(aj->getChainId()) != 0){
            //atom i and atom j are on different chains
            minx=j;
            iMinus6=dist;
          }
          else{
            //atom i and atom j are on the same chain
            diffResId=aj->getResId() - ak->getResId();
            if (diffResId <= -6){
              minx=j;
              iMinus6=dist;
            }
            else{
              if (diffResId == 0 && (Misc::atoi(aj->getICode()) - Misc::atoi(ak->getICode()) <= -6)){
                minx=j;
                iMinus6=dist;
              }
            }
          }
        }
      }

      int k;
      unsigned int q;
      int max;
      max=10;
      if (iatom == 0){
        for (k=0; k< max; k++){
          ai->addData(defVal);
        }
        start=0;
      }
      else{
        start=iatom-1;
      }
      for (q=start; q<= iatom+1; q++){
        if (q < c->getAtmVecSize()){
          ak=c->getAtom(q);
          for (k=static_cast<int>(pinx)-2; k<=static_cast<int>(pinx)+2; k++){
            if (k >= 0 && k < static_cast<int>(natom)){
              aj=pcmol->getAtom(k);
              ai->addData(pcPairDist.at(ak->getAtmInx()).at(aj->getAtmInx()));
            }
            else{
              ai->addData(defVal);
            }
          }
          for (k=static_cast<int>(minx)-2; k<=static_cast<int>(minx)+2; k++){
            if (k >= 0 && k < static_cast<int>(natom)){
              aj=pcmol->getAtom(k);
              ai->addData(pcPairDist.at(ak->getAtmInx()).at(aj->getAtmInx()));
            }
            else{
              ai->addData(defVal);
            } 
          }
        }
        else{
          for (k=0; k< max; k++){
            ai->addData(defVal);
          }
        }
      }

    } //Loop through atoms
  }//Loop through chains

  //Always clear and resize!
  fdataIO.clear();
  fdataIO.resize(camol->getNAtom());

  //Store features
  natom=0; //This is needed since natom is used for other things above
  for (unsigned int ichain=0; ichain < camol->getChnVecSize(); ichain++){
    c=camol->getChain(ichain);
    for (unsigned int iatom=0; iatom < c->getAtmVecSize(); iatom++){
      ai=c->getAtom(iatom);
      fdataIO.at(natom).reserve(3*ai->getDataSize());
  
      //Store S(i)
      for (j=0; j< ai->getDataSize(); j++){
        fdataIO.at(natom).push_back(ai->getDataPoint(j));
      }

      //Store S(i-1)
      if (iatom > 0){
        //Not first atom of chain
        aj=c->getAtom(iatom-1);
        for (j=0; j< ai->getDataSize(); j++){
          if (aj != NULL){
            fdataIO.at(natom).push_back(aj->getDataPoint(j));
          }
          else{
            fdataIO.at(natom).push_back(defVal);
          }
        }
      }
      else{
        //Store S(i-1) which is all default values
        for (j=0; j< ai->getDataSize(); j++){
          fdataIO.at(natom).push_back(defVal);
        }
      }

      //Store S(i+1)
      aj=c->getAtom(iatom+1);
      for (j=0; j< ai->getDataSize(); j++){
        if (aj != NULL){
          fdataIO.at(natom).push_back(aj->getDataPoint(j));
        }
        else{
          fdataIO.at(natom).push_back(defVal);
        }
      }
      natom++;
    } //Loop through atoms
  } //Loop through chains


  if (camol != NULL){
    delete camol;
  }
  if (pcmol != NULL){
    delete pcmol;
  }
}


Eigen::Matrix3d Analyze::gyrationTensor(Molecule* mol){
  Eigen::Matrix3d S;
  Atom* a;
  Coor xyz;
  double dx, dy, dz;

  S=Eigen::Matrix3d::Zero();

  xyz=Analyze::centerOfGeometry(mol);

  //        [SUM(xi-xcog)(xi-xcog) SUM(xi-xcog)*(yi-ycog) SUM(xi-xcog)(zi-zcog)]
  //S=(1/N) [SUM(xi-xcog)(yi-ycog) SUM(yi-ycog)*(yi-ycog) SUM(yi-ycog)(zi-zcog)]
  //        [SUM(xi-xcog)(zi-zcog) SUM(yi-ycog)*(zi-zcog) SUM(zi-zcog)(zi-zcog)]

  for (unsigned int i=0; i< mol->getAtmVecSize(); i++){
    a=mol->getAtom(i);
    dx=a->getX()-xyz.x(); //x(i) - xcog
    dy=a->getY()-xyz.y(); //y(i) - ycog
    dz=a->getZ()-xyz.z(); //z(i) - zcog
    
    S(0,0)+=dx*dx;
    S(0,1)+=dx*dy;
    S(0,2)+=dx*dz;

    S(1,0)+=dy*dx;
    S(1,1)+=dy*dy;
    S(1,2)+=dy*dz;

    S(2,0)+=dz*dx;
    S(2,1)+=dz*dy;
    S(2,2)+=dz*dz;
  }

  S=S/mol->getAtmVecSize();

  return S;
}

double Analyze::quasiharmonicEntropy(Molecule* mol, const Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd>& eigenin, const std::vector<unsigned int> modesin, const double temp){
  unsigned int i, j;
  int nrow;
  Atom* atm;
  double Stot;
  double H, Q, U, F, S, CONST;
  std::vector<double> w; //Frequency
  std::vector<double> wave; //Wavenumbers [1/cm]
  
  nrow=eigenin.eigenvalues().rows();
  Stot=0.0;

  j=0; 
  for (i=0; i< mol->getNAtom(); i++){
    atm=mol->getAtom(i);
    if (atm->getSel() == false){
      continue;
    }
    //Note that the eigenvalues are already in CHARMM units of [amu*Angstroms*Angstroms]
    //So there is no need to multiply it by 10E20/AMU
    if (j < static_cast<unsigned int>(nrow) && std::find(modesin.begin(), modesin.end(),j+1) != modesin.end() && eigenin.eigenvalues()[nrow-(j+1)] > 0){
      w.push_back(sqrt((kB*temp)/(atm->getMass()*eigenin.eigenvalues()[nrow-(j+1)])));
      j++;
    }
    if (j < static_cast<unsigned int>(nrow) && std::find(modesin.begin(), modesin.end(),j+1) != modesin.end() && eigenin.eigenvalues()[nrow-(j+1)] > 0){
      w.push_back(sqrt((kB*temp)/(atm->getMass()*eigenin.eigenvalues()[nrow-(j+1)])));
      j++;
    }
    if (j < static_cast<unsigned int>(nrow) && std::find(modesin.begin(), modesin.end(),j+1) != modesin.end() && eigenin.eigenvalues()[nrow-(j+1)] > 0){
      w.push_back(sqrt((kB*temp)/(atm->getMass()*eigenin.eigenvalues()[nrow-(j+1)])));
      j++;
    }
    //Note that 6 rigid body translation rotations are not explicitly removed!
  }

  wave.resize(w.size());
  CONST=(speedl/PS2SEC)*PLANCK/(KBOLTZ);

  //Andricioaei & Karplus, (2001), J. Chem. Phys. 115:6289
  for (i=0; i< w.size(); i++){
    wave.at(i)=FRQ2INVCM*w.at(i); //Units of 1/cm
    U=CONST*wave.at(i)/temp; //Unitless, no need to multiply wave.at(i) by 1E-2/SPEEDL
    Q=1.0/(1.0-exp(-U));
    F=log(1.0/Q);
    if (U > 40.0){
      //Prevent overflow??
      U=40.0;
    }
    H=U/(exp(U)-1.0);
    S=H-F;
    Stot+=S;
  }
  Stot*=kB; //kcal/mol/K
  
  return Stot;
}

double Analyze::configurationalEntropy(const Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd>& eigenin, const std::vector<unsigned int>& modesin, double cutoffin){
  //Karplus & Kushick, (1981), Macromolecules, 14:325
  int n;
  int nrow;
  double S;
  double det;

  n=eigenin.eigenvalues().rows();
  nrow=eigenin.eigenvalues().rows();
  det=1;
  S=0;

  //for (int i=0; i< nrow; i++){
  for (unsigned int i=0; i< modesin.size(); i++){
    if (modesin.at(i) > static_cast<unsigned int>(nrow)){
      std::cerr << "Warning: Skipping unknown Mode " << modesin.at(i) << std::endl;
    }
    else{
      if (eigenin.eigenvalues()[nrow-modesin.at(i)] > cutoffin){
        det*=eigenin.eigenvalues()[nrow-modesin.at(i)];
      }
    }
  }

  if (det != 0){
    S=(0.5*n*kB)+0.5*kB*log(pow(2*PI, n)*det); //Units are in kcal/mol/K, det is unitless
  }
  else{
    std::cerr << "Warning: Covariance matrix is singular!" << std::endl;
  }

  return S;
}

void AnalyzePcasso::setOutType(PcassoOutEnum pin){
  pout=pin;
}

PcassoOutEnum AnalyzePcasso::getOutType(){
  return pout;
}

