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

#include "Misc.hpp"
#include "Eigen/Eigenvalues"

#include <fstream>

void usage(){
  std::cerr << "Usage:   Covariance [-options] <input(s)>" << std::endl;
  std::cerr << "Options: [-col col1[:col2[:..[:colN]]] || -col col1=colN=incr]]]" << std::endl;
  exit(0);
}

int main (int argc, char **argv){

  int i;
  unsigned int j, k, l;
  std::string currArg;
  std::vector<std::string> ifiles;
  std::ifstream inpFile;
  std::istream* finp;
  std::string line;
  std::vector<double> s;
  std::vector<unsigned int> cols;
  std::vector<unsigned int> range;
  unsigned int incr;
  std::vector<double> avg;
  Eigen::MatrixXd tCovar;
  Eigen::MatrixXd avgCovar;
  unsigned int n;

  ifiles.clear();
  avg.clear();
  n=0;

  for (i=1; i<argc; i++){
    currArg=argv[i];
    if (currArg.compare("-h") == 0 || currArg.compare("-help") == 0){
      usage();
    }
    else if (currArg.compare("-col") == 0 || currArg.compare("-cols") == 0){
      currArg=argv[++i];
      if (currArg.find(":") != std::string::npos || Misc::isdigit(currArg)){
        Misc::splitNum(currArg, ":", cols, false);
      }
      else if (currArg.find("=") != std::string::npos){
        Misc::splitNum(currArg, "=", range, false);
        incr=1;
        if (range.size() >= 2){
          if (range.size() >= 3){
            incr=range.at(2);
          }
        }
        for (j=range.at(0); j<= range.at(1); j=j+incr){
          cols.push_back(j-1);
        }
      }
      else{
        std::cerr << std::endl << "Error: Unrecognized column specification format";
        std::cerr << std::endl << std::endl;
        usage();
      } 
    }
		else if (currArg.compare(0,1,"-") == 0){
      std::cerr << "Warning: Skipping unknown option \"" << currArg << "\"" << std::endl;
    }
    else{
      ifiles.push_back(currArg);
    }
  }

  avg.resize(cols.size());
  for (j=0; j< cols.size(); j++){
    avg.at(j)=0.0;
  }

  tCovar=Eigen::MatrixXd::Zero(cols.size(), cols.size());
  avgCovar=Eigen::MatrixXd::Zero(cols.size(), cols.size());

  //Get average
  for (j=0; j< ifiles.size(); j++){
    inpFile.open(ifiles.at(j).c_str(), std::ios::in);
    finp=&inpFile;

    while(finp->good() && !(finp->eof())){
      getline(*finp, line);
      Misc::splitNum(line, " \t", s, false);
      if (line.length() == 0 || s.size() == 0){
        continue;
      }
      for (k=0; k< cols.size(); k++){
        if (cols.at(k) < s.size()){
          avg.at(k)+=s.at(cols.at(k));
        }
        else{
          std::cerr << "Warning: Missing data in column " << cols.at(k)+1;
          std::cerr << " and line " << n+1 << std::endl;
        }
      }
      n++; 
    }

    if (inpFile.is_open()){
      inpFile.close();
    }
  }

  for (k=0; k< cols.size(); k++){
    avg.at(k)/=n;
  }


  //Get average difference 
  for (j=0; j< ifiles.size(); j++){
    inpFile.open(ifiles.at(j).c_str(), std::ios::in);
    finp=&inpFile;
    while(finp->good() && !(finp->eof())){
      getline(*finp, line);
      Misc::splitNum(line, " \t", s, false);
      if (line.length() == 0 || s.size() == 0){
        continue;
      }
      for (k=0; k< cols.size(); k++){
        tCovar(k,k)=0.0;
        if (cols.at(k) < s.size()){
          tCovar(k,k)=s.at(cols.at(k))-avg.at(k);
        }
      }
      //Generate covariance matrix for all off diagonal elements first
      //Note that the diagonal is NOT zero!
      for (k=0; k< cols.size(); k++){
        for (l=k+1; l< cols.size(); l++){
          tCovar(k,l)=tCovar(k,k)*tCovar(l,l);
          tCovar(l,k)=tCovar(k,l);
        }
        tCovar(k,k)=tCovar(k,k)*tCovar(k,k);
      }
      //Add to average covariance matrix
      avgCovar+=tCovar;
    }

    if (inpFile.is_open()){
      inpFile.close();
    }
  } 

  //Take average
  if (n > 1){
    avgCovar/=(n-1);
  }
  std::cout << avgCovar << std::endl;


  return 0;
}
