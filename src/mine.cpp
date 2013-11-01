//Sean M. Law

#include "MINE/cppmine.h"

#include <iostream>
#include <fstream>
#include <cstdlib>

void usage(){
	std::cerr << std::endl << std::endl;
	std::cerr << "Usage:   mine [-options] <input>" << std::endl;
	std::cerr << "Options: [-x col] [-y col]" << std::endl;
	std::cerr << "         [-mic]" << std::endl;
  exit(0);
}

int main (int argc, char **argv){

  int i;
  std::string currArg;
	std::string inp;

  for (i=1; i<argc; i++){
    currArg=argv[i];
    if (currArg.compare("-h") == 0 || currArg.compare("-help") == 0){
      usage();
    }
    else{
      inp=currArg;
    }
  }

  if (inp.length() == 0){
    std::cerr << std::endl << "Error: Please provide an input file" << std::endl << std::endl;
    usage();
  }

  

  return 0;
}
