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
#include "MINE/cppmine.h"

#include <fstream>
#include <cstdlib>
#include <limits>

void usage(){
  std::cerr << std::endl << std::endl;
  std::cerr << "Usage:   mine [-options] <input>" << std::endl;
  std::cerr << "Options: [-x col] [-y col]" << std::endl;
  std::cerr << "         [-delimiter expression]" << std::endl;
  std::cerr << "         [-start line] [-stop line] [-skip nlines]" << std::endl;
  std::cerr << "         [-max value] [-min value]" << std::endl;
  std::cerr << "         [-mic]" << std::endl;
  exit(0);
}

int main (int argc, char **argv){

  int i;
  std::string currArg;
  std::string inp;
  std::ifstream inpFile;
  std::istream* finp;
  std::string line;
  std::string delimiter;
  std::vector<std::string> s;
  double xval, yval;
  MINE *mine;
  bool mic;
  unsigned int xcol, ycol;
  std::vector<double> xvec, yvec;
  double* x;
  double* y;
  unsigned int nline, start, stop, skip;
  double max, min;

  delimiter=" \t";
  mine=NULL;
  mic=false;
  s.clear();
  xcol=1-1;
  ycol=2-1;
  nline=0;
  start=0;
  stop=std::numeric_limits<unsigned int>::max();
  skip=0;
  min=std::numeric_limits<double>::min();
  max=std::numeric_limits<double>::max();

  for (i=1; i<argc; i++){
    currArg=argv[i];
    if (currArg.compare("-h") == 0 || currArg.compare("-help") == 0){
      usage();
    }
    else if (currArg.compare("-xcol") == 0 || currArg.compare("-x") == 0){
      currArg=argv[++i];
      std::stringstream(currArg) >> xcol;
      xcol--;
    }
    else if (currArg.compare("-ycol") == 0 || currArg.compare("-y") == 0){
      currArg=argv[++i];
      std::stringstream(currArg) >> ycol;
      ycol--;
    }
    else if (currArg.compare("-delimit") == 0 || currArg.compare("-delimiter") == 0){
      currArg=argv[++i];
      delimiter=currArg;
    }
    else if (currArg.compare("-mic") == 0){
      mic=true;
    }
    else if (currArg.compare("-start") == 0){
      currArg=argv[++i];
      std::stringstream(currArg) >> start;
    }
    else if (currArg.compare("-stop") == 0){
      currArg=argv[++i];
      std::stringstream(currArg) >> stop;
    }
    else if (currArg.compare("-skip") == 0){
      currArg=argv[++i];
      std::stringstream(currArg) >> skip;
    }
    else if (currArg.compare("-max") == 0){
      currArg=argv[++i];
      std::stringstream(currArg) >> max;
    }
    else if (currArg.compare("-min") == 0){
      currArg=argv[++i];
      std::stringstream(currArg) >> min;
    }
    else if (currArg.compare(0,1,"-") == 0){
      std::cerr << "Warning: Skipping unknown option \"" << currArg << "\"" << std::endl;
    }
    else{
      inp=currArg;
    }
  }

  if (inp.length() == 0){
    std::cerr << std::endl << "Error: Please provide an input file" << std::endl << std::endl;
    usage();
  }
  else{
    //Read columns into memory
    inpFile.open(inp.c_str(), std::ios::in);
    finp=&inpFile;
    while (finp->good() && !(finp->eof())){
      getline(*finp, line);
      nline++;
      if (line.length() == 0 || nline < start || nline > stop){
        continue;
      }
      else if (skip > 0 && nline % (skip+1) !=0){
        continue;
      }
      else{
        Misc::splitStr(line, delimiter, s, false);
        if (s.size() >= xcol && s.size() >= ycol){
          std::stringstream(s.at(xcol)) >> xval;
          std::stringstream(s.at(ycol)) >> yval;
          if (xval >= min && xval <=max && yval >= min && yval <= max){
            xvec.push_back(xval);
            yvec.push_back(yval);
          }
        }
        else{
          std::cerr << "Warning: Missing columns found and line " << nline << " was ignored" << std::endl;
        }
      }
    }

    if (xvec.size() == yvec.size() && xvec.size() > 0){
      mine = new MINE(0.6, 15);
      x=&xvec[0];
      y=&yvec[0];
      
      mine->compute_score(x, y, static_cast<int>(xvec.size()));
      if (mic == true){
        std::cout << "MIC: " << mine->mic() << std::endl;
      } 
    }
  }

  if (mine != NULL){
    delete mine;
  }
  

  return 0;
}
