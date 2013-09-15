//Sean M. Law

#include "Histogram.hpp"

#include <iostream>
#include <fstream>
#include <cstdlib>

void usage(){
  std::cerr << std::endl;
  std::cerr << "Usage:   histo [-options] <-col colA[:colB[:...[:colN]]]> <inpFile(s)>" << std::endl;
  std::cerr << "Options: [-bins number | -res number]" << std::endl;
  std::cerr << "[-verbose] [-delimiter expression]" << std::endl;
  std::cerr << "[-separate]" << std::endl;
  std::cerr << std::endl;
  exit(0);
}

int main (int argc, char **argv){


  int i;
  unsigned int j, k, m;
  std::vector<Histogram*> histos;
  std::vector<std::string> inps;
  std::string currArg;
  std::vector<double> s;
  std::vector<double> data;
  std::ifstream inpFile;
  std::istream *finp;
  std::vector<std::vector<unsigned int> > cols;
  bool verbose;
  std::string delim;
  std::string line;
  bool separate;

  histos.clear();
  inps.clear();
  finp=NULL;
  cols.clear();
  verbose=false;
  delim=" \t";
  line.clear();
  separate=false;

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
    else if (currArg.compare("-col") == 0 || currArg.compare("-cols") == 0){
      currArg=argv[++i];
      cols.resize(cols.size()+1); //Expand vector
      Misc::splitNum(currArg, ":", cols.at(cols.size()-1), false);
      //Create histogram with one window and nDim = number of columns
      histos.push_back(new Histogram(1, cols.at(cols.size()-1).size()));
    }
    else if (currArg.compare("-verbose") == 0){
      currArg=argv[++i];
      verbose=true;
    }
    else if (currArg.compare("-separate") == 0){
      currArg=argv[++i];
      separate=true;
    }
    else{
      inps.push_back(currArg);
    }
  }


  for (j=0; j< inps.size(); j++){
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
        for (k=0; k< histos.size(); k++){
          data.clear();
          data.resize(cols.at(k).size(),0);
          for (m=0; m< cols.at(k).size(); m++){
            if (cols.at(k).at(m) < s.size()){
              data.at(m)=s.at(cols.at(k).at(m));
            }
            else{
              std::cerr << "Error: Missing data in column ";
              std::cerr << cols.at(k).at(m) << std::endl;
              usage();
            }
          }
          //Append vector of data into to zeroth window
          histos.at(k)->appendData(data, 0);
        }
      }
    }
    if (inps.at(j).compare("-") !=0 && inpFile.is_open()){
      inpFile.close();
    }

    if (separate == true){
      //Output histograms for this file and clear data
      for (k=0; k< histos.size(); k++){
        //histos.setBins(???);
      }
    }
  }

  if (separate == false){
    //Output histograms
    for (k=0; k< histos.size(); k++){
      //histos.setBins(???);
    }
  }

  return 0;
}
