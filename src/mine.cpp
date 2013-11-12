//Sean M. Law

#include "Misc.hpp"
#include "MINE/cppmine.h"

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <vector>

void usage(){
	std::cerr << std::endl << std::endl;
	std::cerr << "Usage:   mine [-options] <input>" << std::endl;
	std::cerr << "Options: [-x col] [-y col]" << std::endl;
	std::cerr << "         [-delimiter expression]" << std::endl;
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
	double tmp;
	MINE *mine;
	bool mic;
	unsigned int xcol, ycol;
	std::vector<double> xvec, yvec;
	double* x;
	double* y;

	delimiter=" \t";
	mine=NULL;
	mic=false;
	s.clear();
	xcol=1-1;
	ycol=2-1;

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
			if (line.length() == 0){
				continue;
			}
			Misc::splitStr(line, delimiter, s, false);
			if (s.size() >= xcol && s.size() >= ycol){
				std::stringstream(s.at(xcol)) >> tmp;
				xvec.push_back(tmp);
				std::stringstream(s.at(ycol)) >> tmp;
				yvec.push_back(tmp);
			}
			else{
				std::cerr << "Warning: Missing columns found and the row was ignored" << std::endl;
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
