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
	std::vector<std::string> tokens;
	std::string delim=":";

/*
	f.resize(3);
	for (unsigned int i=0; i< f.size(); i++){
		f.at(i).resize(3);
		for (unsigned int j=0; j< f.at(i).size(); j++){
			f.at(i).at(j)=(i+1)*(j+1);
//			std::cerr << f.at(i).at(j) << std::endl;
		}
	}
*/

	f.resize(1);
	f.at(0).resize(11);
	f.at(0).at(0)=0.5;
	f.at(0).at(1)=3.5;
	f.at(0).at(2)=4.5;
	f.at(0).at(3)=0.5;
	f.at(0).at(4)=6.0;
	f.at(0).at(5)=5.0;
	f.at(0).at(6)=6.0;
	f.at(0).at(7)=7.0;
	f.at(0).at(8)=10.0;
	f.at(0).at(9)=11.0;
	f.at(0).at(10)=12.0;

//	std::string serialT = "1.0:1 2.0:2 3.0:3 4.0:4 A B 5.0:5 C D 6.0:6 E F 7.0:7 8.0:8 9.0:9 G H 10.0:10 I J 11.0:11 K L";
	std::string serialT = "1.0:1 2.0:2 3.0:3 4.0:4 A B Z 5.0:5 C D 6.0:6 E F 7.0:7 8.0:8 9.0:9 G H 10.0:10 I J 11.0:11 K L";


	Misc::splitStr(Misc::trim(serialT), " \t", tokens, false);
	t->genDTree(tokens, delim);

	std::cerr << t->getDTreeClass(f.at(0)) << std::endl;

	t->delDTree();
/*
	t->addDTree(4.5, 0, "C");
	t->addDTree(5.5, 0, "H");
	t->addDTree(6.5, 0, "E");
	std::cerr << t->getDTreeRoot()->cls << std::endl;
	std::cerr << t->getDTreeRoot()->left->cls << std::endl;
	std::cerr << t->getDTreeRoot()->getDTreeNodeLeft()->cls << std::endl;
	std::cerr << t->getDTreeRoot()->left->left->cls << std::endl;
*/	

  return 0;
}
