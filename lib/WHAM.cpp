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
  metainp=NULL;
  bins.clear();
  tol=1E-5;
  maxIter=1E6;
  T.clear();
  T.push_back(300);
  targetT=0;
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

  line.clear();
  s.clear();

  metaFile.open(fMeta.c_str(), std::ios::in);
  metainp=&metaFile;
  while (metainp->good() && !(metainp->eof())){
    getline(*metainp, line);
    Misc::splitStr(line, " \t", s, false);
    if (s.size() == 3 && (s.at(0) != "!" || s.at(0) != "#")){
      inps.resize(this->getNWindow()+1);
      inps.at(this->getNWindow()).resize(3);
      inps.at(this->getNWindow())=s;
      this->setNWindow(this->getNWindow()+1);
    }
  }

  this->fixTemp(); //Ensure the number of windows and temperatures match, assign targetT
}

void WHAM::processVtot (){
  unsigned int i;
  unsigned int j;

  j=0;
  for(i=0; i< inps.size(); i++){

  }
}

void WHAM::processEtot (){
  unsigned int i;
  unsigned int j;

  j=1;
  for(i=0; i< inps.size(); i++){
    
  }
}

void WHAM::processCoor (){
  unsigned int i;
  unsigned int j;

  j=2;
  for(i=0; i< inps.size(); i++){

  }
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
              denomInv.at(j).at(k)+=nFlast.at(l)*expBVE.at(l).at(j).at(k);
            }
            denomInv.at(j).at(k)=1.0/denomInv.at(j).at(k);
          }
					if (expBVxEx.size() == expBVE.size()){ //WHAM Extrapolation
						FnextInv.at(i)+=expBVxEx.at(i).at(j).at(k)*denomInv.at(j).at(k);
					} 
					else{ //Traditional WHAM
          	FnextInv.at(i)+=expBVE.at(i).at(j).at(k)*denomInv.at(j).at(k);
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
  while (T.size() < this->getNWindow()){
    std::cerr << "Warning: Simulation window " << T.size()+1 << " set to default temperature (";
    std::cerr << targetT << " K)" << std::endl;
    T.push_back(targetT);
  }

  if (T.size() > this->getNWindow()){
    targetT=T.back();
    T.pop_back();
//    std::cerr << T.at(this->getNWindow()-1) << " " << targetT << std::endl;
  }

  if (targetT == 0.0){
    std::cerr << "Warning: Target temperature has been set to its default (300 K)" << std::endl;
    targetT=300;
  }
  
  if (T.size() > this->getNWindow()){
    std::cerr << "Warning: Extra temperature(s) provided were ignored" << std::endl;
    T.resize(this->getNWindow());
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
  T.clear();
  Misc::splitNum(tin, ":", T);
  if (T.size() <= 0){
    std::cerr << std::endl << "Error: Unrecognized temperature format " << tin;
    std::cerr << std::endl << std::endl;
    return true;
  }
  else{
    return false;
  }
}

void WHAM::setTemp(const std::vector<double> &tin){
  T.clear();
  T=tin;
}

bool WHAM::setTempRange(const std::string &tin){
  std::vector<double> s;
  unsigned int i;

  T.clear();
  Misc::splitNum(tin, "=", s);

  if (s.size() >= 3){
    for (i=0; i<= static_cast<unsigned int>((s.at(1)-s.at(0))/s.at(2)); i++){
      T.push_back(s.at(0)+i*s.at(2));
    }
    if (s.size() >=4){
      T.push_back(s.at(3)); //Target temp
    }
    return false;
  }
  else if (s.size() == 2){
    for (i=0; i<= static_cast<unsigned int>((s.at(1)-s.at(0))); i++){
      T.push_back(s.at(0)+i);
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
  return T.size();
}

double WHAM::getTemp(const int &element){
  return T.at(element);
}

std::string WHAM::getCmd(){
  return cmd;
}

unsigned int WHAM::getNWindow(){
  return nWindow;
}
