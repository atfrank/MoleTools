//Sean M. Law

#include "Histogram.hpp"

#include <iostream>
#include <fstream>
#include <cstdlib>

void usage(){
  std::cerr << "Usage:   histo [-options] <inpFile(s)>" << std::endl;
  std::cerr << "Options: [-bins number | -res number]" << std::endl;
  std::cerr << "[-col colA[:colB[:...[:colN]]]]" << std::endl;
  std::cerr << "[-verbose] [-delimiter expression]" << std::endl;
  std::cerr << std::endl;
  exit(0);
}

int main (int argc, char **argv){


  int i;
  std::vector<Histogram*> histos;
  std::vector<std::string> inps;
  std::string currArg;
  std::vector<double> s;
  std::ifstream inpFile;
  std::istream *finp;
  std::vector<std::vector<unsigned int> > col;
  bool verbose;
  std::string delim;
  std::string line;

  histos.clear();
  inps.clear();
  finp=NULL;
  col.clear();
  verbose=false;
  delim=" \t";
  line.clear();

  for (i=1; i<argc; i++){
    currArg=argv[i];
    if (currArg.compare("-h") == 0 || currArg.compare("-help") == 0){
      usage();
    }
    else if (currArg.compare("-bins") == 0){
      currArg=argv[++i];
    }
    else if (currArg.compare("-res") == 0){
      //Resolution/Bin spacing

    }
    else if (currArg.compare("-delimiter") == 0 || currArg.compare("-delimit") == 0){
      currArg=argv[++i];
      delim=currArg;
    }
    else{
      inps.push_back(currArg);
    }
  }

  for (unsigned int j=0; j< inps.size(); j++){
    if (verbose == true){
      if (inps.at(j).compare("-") == 0){
        std::cerr << "Processing data from STDIN" << std::endl;
      }
      else{
        std::cerr << "Processing file \"" << inps.at(j) << "\"..." << std::endl;
      }
    }
    if (inps.at(j).compare("-") == 0){
      finp=&std::cin;
    }
    else{
      inpFile.open(inps.at(j).c_str(), std::ios::in);
      finp=&inpFile;
    }
    while (finp->good() && !(finp->eof())){
      getline(*finp, line);
      if (line.length() > 0){
        Misc::splitNum(line, delim, s, false);
        
      }
    }
    if (inps.at(j).compare("-") !=0 && inpFile.is_open()){
      inpFile.close();
    }
  }

  return 0;
}
