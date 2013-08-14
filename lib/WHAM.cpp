// Sean M. Law

#include "WHAM.hpp"

WHAM::WHAM (){
  cmd.clear();
  Fguess.clear();
  F.clear();
  expBVE.clear();
	expBVxEx.clear();
  nWindow=0;
  fMeta.clear();
  bins.clear();
  tol=1E-5;
  maxIter=1E6;
  B.clear();
  B.push_back(1.0/(kB*300));
  B0=0;
  factor=1.0;
  factorFlag=false;
}

void WHAM::appendCmd(const std::string &str){
  cmd+=str;
}

void WHAM::genWHAMInput(){

}

void WHAM::readMetadata(){
  std::string line;
  std::vector<std::string> s;
  std::ifstream metaFile;
  std::istream *metainp;


  line.clear();
  s.clear();
  metainp=NULL;

  metaFile.open(fMeta.c_str(), std::ios::in);
  metainp=&metaFile;
  while (metainp->good() && !(metainp->eof())){
    getline(*metainp, line);
    Misc::splitStr(line, " \t", s, false);
    if (s.size() == 3 && (s.at(0) != "!" || s.at(0) != "#")){
      inps.resize(this->getNWindow()+1);
      inps.at(this->getNWindow())=s;
      this->setNWindow(this->getNWindow()+1);
    }
  }

  if (metaFile.is_open()){
    metaFile.close();
  }

  this->fixTemp(); //Ensure the number of windows and temperatures match, assign B0 
}

void WHAM::processEnergies(){
  unsigned int i; //Simulation environment
  unsigned int j; //Simulation window
  unsigned int k; //Datapoint in simulation window j
  std::string vline; //Biasing potential
  std::string eline; //Total potential energy
  std::vector<double> v;
  std::vector<double> e;
  std::ifstream fvin;
  std::ifstream fein;
  std::istream *fvinp;
  std::istream *feinp;
  unsigned int lastVSize;

  fvinp=NULL;
  feinp=NULL;
  lastVSize=0;

  expBVE.resize(inps.size());

  for(j=0; j< inps.size(); j++){
    fvin.open(inps.at(j).at(0).c_str(), std::ios::in);
    fein.open(inps.at(j).at(1).c_str(), std::ios::in);
    fvinp=&fvin;
    feinp=&fein;
    k=0; //Datapoint in simulation window j

    if (fvinp->good() && feinp->good()){
      //Read both files
      while (fvinp->good() && !(fvinp->eof()) && feinp->good() && !(feinp->eof())){
        getline(*fvinp, vline);
        getline(*feinp, eline);
        Misc::splitNum(vline, " \t", v, false);
        Misc::splitNum(eline, " \t", e, false);
        if (vline.length() ==0){
          getline(*fvinp, vline);
          if (fvinp->eof()){
            continue;
          }
          else{
            std::cerr << "Error: file contains too many lines" << std::endl;
          }
        }
        if (eline.length() == 0){
          getline(*feinp, eline);
          if (feinp->eof()){
            continue;
          }
          else{
            std::cerr << "Error: file contains too many lines" << std::endl;
          }
        }
        expBVE.at(j).resize(k+1);
        if (k == 0){
          lastVSize=v.size();
        }
        if (e.size() == 1 && this->getNWindow() == this->getTempSize()){
          if (v.size() == this->getNWindow() && v.size() == lastVSize){
            //Traditional WHAM
            //Consider using std::transform instead
            for (i=0; i< v.size(); i++){
              v.at(i)=exp(-B.at(i)*v.at(i));
            }
            expBVE.at(j).at(k)=v;
            for (i=0; i< e.size(); i++){
              expBVE.at(j).at(k).at(i)*=exp(-(B.at(i)-B0)*e.at(i));
            }
          }
          else if (v.size() == 2*this->getNWindow() && v.size() == lastVSize){
          //WHAM Extrapolation

          }
          else{
            std::cerr << "Warning: File \"???\" line "<< k+1;
            std::cerr << " contains the wrong number of columns" << std::endl;
            continue;
          }
        }
        k++;
      }
    }
    else if (fvinp->good()){
      //Read biasing potential only
      while (fvinp->good() && !(fvinp->eof())){
        getline(*fvinp, vline);
        Misc::splitNum(vline, " \t", v, false);
        if (k == 0){
          lastVSize=v.size();
        }
      }
    }
    else if (feinp->good()){
      //Read total potential only
      while (feinp->good() && !(feinp->eof())){
        getline(*feinp, eline);
        Misc::splitNum(eline, " \t", e, false);
        
      }
    }
    else{
      std::cerr << "Warning: Files \"" << inps.at(i).at(0) << "\" and \"";
      std::cerr << inps.at(i).at(1) << "\" could not be read" << std::endl;
    }

    if (fvin.is_open()){
      fvin.close();
    }
    if (fein.is_open()){
      fein.close();
    }
  }

}


void WHAM::processCoor (){
  
}

bool WHAM::iterateWHAM (){
  unsigned int i,j,k,l;
  unsigned int niter;
  std::vector<double> nFlast;
  std::vector<double> FnextInv;
  std::vector< std::vector<double> > denomInv;
  bool nextIter;
  double dExpF; //exp(F)-exp(Flast)
 
  // WHAM Formalism (Adapted from Michael Andrec)
  //
  // 1.0/F(i) = 
  //
  //                                             exp(-B(i)V(i,jk))*exp[-(B(i)-B(0))E(i,jk)]
  // Sum(j=1,..,S) Sum(k=1,...,N(j)) --------------------------------------------------------------------
  //                                  Sum(l=1,..,S) N(l)*F(l)*exp(-B(l)V(l,jk))*exp[-(B(l)-B(0))E(l,jk)]
  //  

  // Comments
  // expBVE.size() = Number of simulations/windows, S = nWindow
	// expBVE.at(i).size() = Number of datapoints, N(i) in the ith simulation/window
  // expBVE.at(i).at(j).at(k) = exp(-B(i)V(i,jk))*exp[-(B(i)-B(0))E(i,jk)]
 
  if (expBVE.size() != this->getNWindow() || (expBVxEx.size() > 0 && expBVxEx.size() != this->getNWindow())){
    std::cerr << std::endl;
    std::cerr << "Error: The number of simulation windows do not match!";
    std::cerr << std::endl << std::endl;
    return true;
  }

  nFlast.resize(this->getNWindow()); //n(i)*F(i)
  FnextInv.resize(this->getNWindow());
  denomInv.resize(this->getNWindow());

  for (i=0; i< this->getNWindow(); i++){
    //Initialize F
    nFlast.at(i)=expBVE.at(i).size();
    if (Fguess.size() > 0 && Fguess.size() == this->getNWindow()){
      nFlast.at(i)*=Fguess.at(i);
    }
    FnextInv.at(i)=0;
    
    denomInv.at(i).resize(expBVE.at(i).size());
  }

  //WHAM Iterations
  for (niter=0; niter< maxIter; niter++){
    for (i=0; i< this->getNWindow(); i++){ //For each F value I
			FnextInv.at(i)=0.0;
      for (j=0; j< this->getNWindow(); j++){ //For each simulation J
        for (k=0; k< expBVE.at(j).size(); k++){ //Foreach datapoint K in simulation J
          if (i==0){ //Calculate (redundant) denominator once for each iteration
						denomInv.at(j).at(k)=0.0;
            for (l=0; l< this->getNWindow(); l++){ //Foreach simulation environment L
              //Calculate denom
              denomInv.at(j).at(k)+=nFlast.at(l)*expBVE.at(j).at(k).at(l);
            }
            denomInv.at(j).at(k)=1.0/denomInv.at(j).at(k);
          }
					if (expBVxEx.size() == expBVE.size()){ //WHAM Extrapolation
						FnextInv.at(i)+=expBVxEx.at(j).at(k).at(i)*denomInv.at(j).at(k);
					} 
					else{ //Traditional WHAM
          	FnextInv.at(i)+=expBVE.at(j).at(k).at(i)*denomInv.at(j).at(k);
					}
        }
      }
      if (factorFlag == true){ //For use with Molecular Transfer Model (MTM)
        FnextInv.at(i)*=factor;
      }
    }

    //Check tolerance (note that tolerance is in f but F is in exp(f))
    nextIter=false;
    for (i=0; i< this->getNWindow(); i++){
      dExpF=fabs(exp(1.0/FnextInv.at(i)) - exp(nFlast.at(i)/expBVE.at(i).size()));
      if (dExpF >= tol){
        nextIter=true;
      }
    }

    //Update F values for next iteration
    for (i=0; i< this->getNWindow(); i++){
      nFlast.at(i)=expBVE.at(i).size()*(1.0/FnextInv.at(i)); //n(i)*F(i)
    }
  }
  
  return false;
}

void WHAM::fixTemp(){
  while (B.size() < this->getNWindow()){
    std::cerr << "Warning: Simulation window " << B.size()+1 << " set to default temperature (";
    std::cerr << B0 << " K)" << std::endl;
    B.push_back(B0);
  }

  if (B.size() > this->getNWindow()){
    B0=B.back();
    B.pop_back();
//    std::cerr << B.at(this->getNWindow()-1) << " " << B0 << std::endl;
  }

  if (B0 == 0.0){
    std::cerr << "Warning: Target temperature has been set to its default (300 K)" << std::endl;
    B0=1.0/(kB*300);
  }
  
  if (B.size() > this->getNWindow()){
    std::cerr << "Warning: Extra temperature(s) provided were ignored" << std::endl;
    B.resize(this->getNWindow());
  }
 
}

void WHAM::setMeta(const std::string &metain){
  fMeta=metain;
}

void WHAM::setBins(const std::string &binsin){
  Misc::splitNum(binsin, ":", bins);
}

void WHAM::setBins(const std::vector<unsigned int> &binsin){
  bins=binsin;
}

void WHAM::setBins(const std::vector<int> &binsin){
  unsigned int i;
  bins.resize(binsin.size());
  for (i=0; i< binsin.size(); i++){
    bins.at(i)= static_cast<unsigned int>(binsin.at(i));
  }
}

void WHAM::setTol(const double &tolin){
  tol=tolin;
}

void WHAM::setMaxIter(const unsigned int &iterin){
  maxIter=iterin;
}

bool WHAM::setTemp(const std::string &tin){
  B.clear();
  Misc::splitNum(tin, ":", B);
  if (B.size() <= 0){
    std::cerr << std::endl << "Error: Unrecognized temperature format " << tin;
    std::cerr << std::endl << std::endl;
    return true;
  }
  else{
    for (unsigned int i=0; i< B.size(); i++){
      B.at(i)=1.0/(kB*B.at(i));
    }
    return false;
  }
}

void WHAM::setTemp(const std::vector<double> &tin){
  B.clear();
  B=tin;
  for (unsigned int i=0; i< B.size(); i++){
    B.at(i)=1.0/(kB*B.at(i));
  }
}

bool WHAM::setTempRange(const std::string &tin){
  std::vector<double> s;
  unsigned int i;

  B.clear();
  Misc::splitNum(tin, "=", s);

  if (s.size() >= 3){
    for (i=0; i<= static_cast<unsigned int>((s.at(1)-s.at(0))/s.at(2)); i++){
      B.push_back(s.at(0)+i*s.at(2));
    }
    if (s.size() >=4){
      B.push_back(s.at(3)); //Target temp
    }
    for (i=0; i< B.size(); i++){
      B.at(i)=1.0/(kB*B.at(i));
    }
    return false;
  }
  else if (s.size() == 2){
    for (i=0; i<= static_cast<unsigned int>((s.at(1)-s.at(0))); i++){
      B.push_back(s.at(0)+i);
    }
    for (i=0; i< B.size(); i++){
      B.at(i)=1.0/(kB*B.at(i));
    }
    return false;
  }
  else{
    std::cerr << std::endl << "Error: Unrecognized temperature format " << tin;
    std::cerr << std::endl << std::endl;
    return true;
  }
}

void WHAM::setFactor(const double &factorin){
  factorFlag=true;
  factor=factorin;
}

void WHAM::setNWindow(const unsigned int &nwin){
  nWindow=nwin;
}

void WHAM::setNWindow(const int &nwin){
  nWindow=static_cast<unsigned int>(nwin);
}

std::string WHAM::getMeta(){
  return fMeta;
}

unsigned int WHAM::getTempSize(){
  return B.size();
}

double WHAM::getTemp(const int &element){
  return 1.0/(kB*B.at(element));
}

std::string WHAM::getCmd(){
  return cmd;
}

unsigned int WHAM::getNWindow(){
  return nWindow;
}
