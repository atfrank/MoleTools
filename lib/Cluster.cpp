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

#include "Cluster.hpp"

#include <iostream>
#include <fstream>

Cluster::Cluster(){

}

Cluster::Cluster(unsigned int Nin){
  this->setN(Nin);
}

void Cluster::readDMatrix(std::string finp){
  std::ifstream inpFile;
  std::istream* inp;
  std::string line;
  unsigned int nline;
  std::vector<double> s;

  inp=NULL;
  nline=0;

  if (finp.length() == 0){
    std::cerr << std::endl << "Error: Please provide an input file of pairwise distances" << std::endl << std::endl;
    return;
  }

  if (finp.compare("-") == 0){
    inp=&std::cin;
  }
  else{
    inpFile.open(finp.c_str());
    inp=&inpFile;
  }
  while (inp->good() && !(inp->eof())){
    getline(*inp, line);
    if (line.length() > 0){
      nline++;
      continue;
    }
  }

}

void Cluster::setN(unsigned int Nin){
  N=Nin;
  cluster.resize(Nin);
  dmatrix.resize(Nin);
  for (unsigned int i=0; i< Nin; i++){
    dmatrix.at(i).resize(Nin);
  }
}

unsigned int Cluster::getN(){
  return N;
}

void Cluster::preCluster(){

}

void Cluster::runCluster(){
  //Do Nothing
}

void DPCluster::preCluster(){
  
}

void DPCluster::runCluster(){

}

void DPCluster::setDc(const double& dcin){
  dc=dcin;
}
    
double DPCluster::getDc(){
  return dc;
}

void DPCluster::setDmax(const double& dmaxin){
  dmax=dmaxin;
}

double DPCluster::getDmax(){
  return dmax;
}

void DPCluster::setMaxDelta(const double& maxdeltain){
  maxDelta=maxdeltain;
}

double DPCluster::getMaxDelta(){
  return maxDelta;
}

void DPCluster::setRho(){
  this->rho.resize(this->getN());

}

void DPCluster::setDelta(){
  this->delta.resize(this->getN());
  this->setMaxDelta(0.0);
}

void KCluster::runCluster(){

}
