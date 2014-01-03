//Sean M. Law

#include "DTree.hpp"

#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <cmath>
#include <vector>

void usage(){
  std::cerr << std::endl << std::endl;
  std::cerr << "Usage:   unitTest [-options] <file>" << std::endl;
  std::cerr << "Options: [-format type]" << std::endl;
  std::cerr << "         [-top file] [-prm file]" << std::endl;
  exit(0);
}

int main (int argc, char **argv){


  int i;
  std::string currArg;
  std::vector<std::string> s;

  for (i=1; i<argc; i++){
    currArg=argv[i];
    if (currArg.compare("-h") == 0 || currArg.compare("-help") == 0){
      usage();
    }
    else{

    }
  }

	DTree* root=new DTree; //Keep root for this tree
	DTree* t=root;

	t->addDTree(4.5);
	t->searchDTree(4.5)->SSE="C";
	std::cerr << t->searchDTree(4.5)->SSE;
	std::cerr << t->searchDTree(4.5)->left;
	t->searchDTree(4.5)->col=10;
	std::cerr << t->searchDTree(4.5)->col;
	std::cerr << std::endl;

  return 0;
}
