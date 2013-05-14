//Sean M. Law

//Contact Appearance Order (CAO)

#include "Misc.hpp"

#include <iostream>
#include <fstream>
#include <cstdlib>

#define MAXINPFILES 4096
#define MAXLINESIZE 4096

void usage(){
	std::cerr << std::endl << std::endl;
	std::cerr << "Usage:   contactAppOrder [-options] <file(s)>" << std::endl;
	std::cerr << "Options: [-skip lines]" << std::endl;
  exit(0);
}

int main (int argc, char **argv){


  int i,l;
	unsigned int j, k;
  std::string currArg;
	std::vector<std::string> ifiles;
	std::ifstream inpFile;
	std::istream* inp;
	std::string line;
	int skip;
	std::vector<int> s;
	std::vector<int> order;
	std::vector<int> rank;
	unsigned int time;
	int nline;
	int nevents;
	bool trackFlag;

	ifiles.clear();
	skip=0;
	nline=0;
	order.clear();
	time=0;
	nevents=0;
	trackFlag=false;

  for (i=1; i<argc; i++){
    currArg=argv[i];
    if (currArg == "-h" || currArg == "-help"){
      usage();
    }
    else if (currArg == "-skip"){
      currArg=argv[++i];
			std::stringstream(currArg) >> skip;
			if (skip > 0){
				skip--;
			}
    }
    else{
			ifiles.push_back(currArg);
    }
  }

	//Loop through each file
	for (j=0; j< ifiles.size(); j++){
		if (ifiles.at(j) == "-"){
    	inp=&std::cin;
  	}
  	else{
    	inpFile.open((ifiles.at(j)).c_str());
    	inp=&inpFile;
  	}
		while (inp->good() && !(inp->eof())){
			getline(*inp, line);
			if (line.size() == 0){
				continue;
			}
			if (nline != skip){
				nline++;
				continue;
			}
			else{
				//Assess snapshot
				Misc::splitNum(line, " \t", s, false); //Split on one or more consecutive whitespaces
				if (order.size() == 0){
					order.resize((s.size()-3)/2);
				}
				if (rank.size() == 0){
					rank.resize((s.size()-3)/2);
				}
				if (s.at(0) == 0 && trackFlag == false){
					time=0;
					std::fill(order.begin(), order.end(), 0);
					trackFlag=true;
				}
				else if (trackFlag == true){
					time++;
					l=0;
					for (k=3; k< s.size(); k=k+2){
						if (s.at(k) == 1 && order.at(l) == 0){
							order.at(l)=time; //Only record time if the contact was broken
						}
						l++;
					}
				}
				else{
					//Clean up
					if (s.at(0) == s.at(2)){ //All contacts are formed
						trackFlag=false;
						nevents++;
					}
				}
				nline=0;
			}
		}
	}

		

  return 0;
}
