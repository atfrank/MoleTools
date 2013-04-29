//Sean M. Law

#include "Misc.hpp"
#include "Molecule.hpp"

#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <cmath>

#define MAXINPFILES 4096
#define MAXLINESIZE 4096

void usage(){
  std::cerr << std::endl << std::endl;
  std::cerr << "Usage:   unitTest [-options] <file>" << std::endl;
  std::cerr << "Options: " << std::endl;
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

  ifile.clear();

  for (i=1; i<argc; i++){
    currArg=argv[i];
    if (currArg == "-h" || currArg == "-help"){
      usage();
    }
    else{
      ifile=currArg;
    }
  }

  if (ifile.length() == 0){
    std::cerr << std::endl << "Error: Please provide an input file" << std::endl << std::endl;
    usage();
  }

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
    s=Misc::split(line, " \t", false); //Split on one or more consecutive whitespace
    for (j=0; j< s.size(); j++){
  //    std::cerr << s.at(j) << ":";
    }
  //  std::cerr << std::endl;
  }

  if (ifile != "-"){
    inpFile.close();
  }

	Molecule* mol;
	mol=new Molecule;
	mol->readPDB(ifile);
	mol->select(".CA");

  return 0;
}
