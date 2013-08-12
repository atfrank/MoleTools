// Sean M. Law

#include "WHAM.hpp"

WHAM::WHAM (){
  Fguess.clear();
  F.clear();
  E.clear();
	Ex.clear();
  fMeta.clear();
  bins.clear();
  tol=1E-5;
  maxIter=1E6;
}

void WHAM::genWHAMInput(){

}

void WHAM::readWHAMInput(){

}

void WHAM::iterateWHAM (){
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
  // E.size() = Number of simulations/windows, S
	// E.at(i).size() = Number of datapoints, N(i) in the ith simulation/window
  // E.at(i).at(j).at(k) = exp(-B(i)V(i,jk))*exp[-(B(i)-B(0))E(i,jk)]
  

  nFlast.resize(E.size()); //n(i)*F(i)
  FnextInv.resize(E.size());
  denomInv.resize(E.size());

  for (i=0; i< E.size(); i++){
    //Initialize F
    nFlast.at(i)=E.at(i).size();
    if (Fguess.size() > 0 && Fguess.size() == E.size()){
      nFlast.at(i)*=Fguess.at(i);
    }
    FnextInv.at(i)=0;
    
    denomInv.at(i).resize(E.at(i).size());
  }

  //WHAM Iterations
  for (niter=0; niter< maxIter; niter++){
    for (i=0; i< E.size(); i++){ //For each F value I
			FnextInv.at(i)=0.0;
      for (j=0; j< E.size(); j++){ //For each simulation J
        for (k=0; k< E.at(j).size(); k++){ //Foreach datapoint K in simulation J
          if (i==0){ //Calculate (redundant) denominator once for each iteration
						denomInv.at(j).at(k)=0.0;
            for (l=0; l< E.size(); l++){ //Foreach simulation environment L
              //Calculate denom
              denomInv.at(j).at(k)+=nFlast.at(l)*E.at(l).at(j).at(k);
            }
            denomInv.at(j).at(k)=1.0/denomInv.at(j).at(k);
          }
					if (Ex.size() == E.size()){ //WHAM Extrapolation
						FnextInv.at(i)+=Ex.at(i).at(j).at(k)*denomInv.at(j).at(k);
					} 
					else{ //Traditional WHAM
          	FnextInv.at(i)+=E.at(i).at(j).at(k)*denomInv.at(j).at(k);
					}
        }
      }
    }

    //Check tolerance (note that tolerance is in f but F is in exp(f))
    nextIter=false;
    for (i=0; i< E.size(); i++){
      dExpF=fabs(exp(1.0/FnextInv.at(i)) - exp(nFlast.at(i)/E.at(i).size()));
      if (dExpF >= tol){
        nextIter=true;
      }
    }

    //Update F values for next iteration
    for (i=0; i< E.size(); i++){
      nFlast.at(i)=E.at(i).size()*(1.0/FnextInv.at(i)); //n(i)*F(i)
    }
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
