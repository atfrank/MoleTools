// Sean M. Law

#include "WHAM.hpp"

WHAM::WHAM (){
  nSims=0;
  nData.clear();
  Fguess.clear();
  energy.clear();
  bins=10;
  tol=1E-5;
  maxIter=100000;
  extrapFlag=false;
  guessFlag=false;
}

void WHAM::genWHAMInput(){

}

void WHAM::readWHAMInput(){

}

void WHAM::iterateWHAM (){
  unsigned int i,j,k,l;
  unsigned int niter;
  std::vector<double> Flast;
  std::vector<double> FnextInv;
  std::vector< std::vector<double> > denom;
  

  Flast.resize(nSims);
  FnextInv.resize(nSims);
  denom.resize(nSims);
  for (i=0; i< nSims; i++){
    //Initialize F
    if (guessFlag == true){
      Flast.at(i)=Fguess.at(i);
    }
    else{
      Flast.at(i)=1;
    }
    FnextInv.at(i)=0;
    
    denom.at(i).resize(nData.at(i));
  }

  for (niter=0; niter< maxIter; niter++){
    //Clear denominators before each iteration
    for (j=0; j< nSims; j++){
      for (k=0; k< nData.at(j); k++){
        
      }
    }

    for (i=0; i< nSims; i++){ //For each F value I
      for (j=0; j< nSims; j++){ //For each simulation J
        for (k=0; k< nData.at(j); k++){ //Foreach datapoint K in simulation J
          if (i==0){ //Calculate (redundant) denominator once for each iteration
            for (l=0; l< nSims; l++){ //Foreach simulation environment L
              if (extrapFlag == true){ //WHAM Extrapolation
              
              }
              else{ //Traditional WHAM
              
              }
            }
          }
          FnextInv.at(i)+=energy.at(i).at(j).at(k);
        }
      }
    }

    //Update F values

    for (i=0; i< nSims; i++){
      Flast.at(i)=1.0/FnextInv.at(i);
      //Note that tolerance is in f but F is in exp(f)!!
    }
  }

}
