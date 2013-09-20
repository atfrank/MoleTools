//Sean M. Law

#include "Misc.hpp"
#include "Molecule.hpp"
#include "Analyze.hpp"
#include "Prmtop.hpp"

#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <cmath>

void usage(){
  std::cerr << std::endl << std::endl;
  std::cerr << "Usage:   unitTest [-options] <file>" << std::endl;
  std::cerr << "Options: [-format type]" << std::endl;
  std::cerr << "         [-top file] [-prm file]" << std::endl;
  exit(0);
}

int main (int argc, char **argv){


  int i;
  unsigned int j;
//  int xcol;
//  int ycol;
//  int skip;
  std::string ifile;
  std::string currArg;
  std::ifstream inpFile;
  std::istream* inp;
  std::string line;
  std::vector<std::string> s;
  std::vector<double> d;
	std::string format;
//  std::string prm;
  std::string top;

  ifile.clear();
	format.clear();
  inp=NULL;
  j=0;

  for (i=1; i<argc; i++){
    currArg=argv[i];
    if (currArg.compare("-h") == 0 || currArg.compare("-help") == 0){
      usage();
    }
		else if (currArg.compare("-format") == 0){
			currArg=argv[++i];
			Misc::toupper(currArg);
			format=currArg;
		}
    else if (currArg.compare("-top") == 0){
      currArg=argv[++i];
      top=currArg;
    }
    else{
      ifile=currArg;
    }
  }

  if (ifile.length() == 0){
    std::cerr << std::endl << "Error: Please provide an input file" << std::endl << std::endl;
    usage();
  }
/*
  if (ifile == "-"){
    inp=&std::cin; 
  }
  else{
    inpFile.open((ifile).c_str());
    inp=&inpFile;
  }

  while (inp->good() && !(inp->eof())){
    getline(*inp, line);
    if (line.size() == 0){
      continue;
    }
    Misc::splitStr(line, " \t", s, false); //Split on one or more consecutive whitespace
    Misc::splitNum(line, " \t", d, false);
//    Misc::splitNum(line, ":", d,true);  
    for (j=0; j< s.size(); j++){
  //    std::cerr << s.at(j) << ":";
    }
  //  std::cerr << std::endl;
  }

  if (ifile != "-"){
    inpFile.close();
  }
*/

	Molecule* mol;
	mol=Molecule::readPDB(ifile, format);
  Prmtop* t;
  t=new Prmtop;
  t->readTopology(top);
  mol->setMass(t);
  mol->setCharge(t);
  for (unsigned int k=0; k< mol->getAtmVecSize(); k++){
    Atom* a;
    a=mol->getAtom(k);
//    std::cerr << a->getSummary() << " " << a->getMass() << " " << a->getCharge() << std::endl;
  }
	/*
	mol->select("B:.CA");
	mol->storeSel();
	mol->selAll();
	mol->recallSel();
	*/
	//mol->writePDB();

/*
	Analyze* a;
	AnalyzeDistance* dist=new AnalyzeDistance;
	AnalyzeDihedral* phi=new AnalyzeDihedral;
	phi->preAnalysis();
	dist->preAnalysis();
	a=phi;
	a->preAnalysis();
*/
  return 0;
}
