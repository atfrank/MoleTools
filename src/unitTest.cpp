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
	std::vector<std::vector <double> > f;

	f.resize(3);
	for (unsigned int i=0; i< f.size(); i++){
		f.at(i).resize(3);
		for (unsigned int j=0; j< f.at(i).size(); j++){
			f.at(i).at(j)=(i+1)*(j+1);
//			std::cerr << f.at(i).at(j) << std::endl;
		}
	}

/*
	t->addDTree(4.5);
	t->searchDTree(4.5)->SSE="C";
	std::cerr << t->searchDTree(4.5)->SSE << std::endl;
	std::cerr << t->searchDTree(4.5)->left << std::endl;
	t->searchDTree(4.5)->inx=10;
	std::cerr << t->searchDTree(4.5)->inx << std::endl;
	t->addDTree(5.5);
	t->searchDTree(5.5)->inx=20;
	std::cerr << t->searchDTree(4.5)->inx << std::endl;
	std::cerr << t->searchDTree(5.5)->inx << std::endl;
	std::cerr << t->searchDTree(4.5)->right->inx << std::endl;
*/

	t->addDTree(4.5, 0, "C");
	t->addDTree(5.5, 0, "H");
	t->addDTree(6.5, 0, "E");
	std::cerr << t->getDTreeRoot()->cls << std::endl;
	std::cerr << t->getDTreeRoot()->left->cls << std::endl;
	std::cerr << t->getDTreeRoot()->getDTreeNodeLeft()->cls << std::endl;
	std::cerr << t->getDTreeRoot()->left->left->cls << std::endl;
	

  return 0;
}
