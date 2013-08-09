// Sean M. Law

#include "WHAM.hpp"

WHAM::WHAM (){
  Fguess.clear();
  E.clear();
	Ex.clear();
  bins=10;
  tol=1E-5;
  maxIter=100000;
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
  
	//Notes
	// E.size() = Number of simulations/windows
	// E.at(i).size() = Number of datapoints in each simulation/window

  Flast.resize(E.size());
  FnextInv.resize(E.size());
  denom.resize(E.size());

  for (i=0; i< E.size(); i++){
    //Initialize F
    if (Fguess.size() > 0 && Fguess.size() == E.size()){
      Flast.at(i)=Fguess.at(i);
    }
    else{
      Flast.at(i)=1;
    }
    FnextInv.at(i)=0;
    
    denom.at(i).resize(E.at(i).size());
  }

  for (niter=0; niter< maxIter; niter++){
    for (i=0; i< E.size(); i++){ //For each F value I
			FnextInv.at(i)=0.0;
      for (j=0; j< E.size(); j++){ //For each simulation J
        for (k=0; k< E.at(j).size(); k++){ //Foreach datapoint K in simulation J
          if (i==0){ //Calculate (redundant) denominator once for each iteration
						denom.at(j).at(k)=0.0;
            for (l=0; l< E.size(); l++){ //Foreach simulation environment L
              //Calculate denom

            }
          }
					if (Ex.size() == E.size()){ //WHAM Extrapolation
						FnextInv.at(i)+=Ex.at(i).at(j).at(k)/denom.at(j).at(k);
					} 
					else{ //Traditional WHAM
          	FnextInv.at(i)+=E.at(i).at(j).at(k)/denom.at(j).at(k);
					}
        }
      }
    }

    //Update F values

    for (i=0; i< E.size(); i++){
      Flast.at(i)=1.0/FnextInv.at(i);
      //Note that tolerance is in f but F is in exp(f)!!
    }
  }

}
