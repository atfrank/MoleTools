//Sean M. Law

#include "Histogram.hpp"
#include "Misc.hpp"

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <ctime>

void usage(){
  std::cerr << std::endl;
  std::cerr << "Usage:   histo [-options] <-col colA[:colB[:...[:colN]]]> <inpFile(s)>" << std::endl;
  std::cerr << "Options: [-bins bA[:bB[:...[bN]]] | -res rA:[:rB[:...[rN]]]]" << std::endl;
  std::cerr << "         [-verbose] [-delimiter expression] [-separate]" << std::endl;
  std::cerr << "         [-prob | -density] [-energy temp]" << std::endl;
  std::cerr << "         [-start line] [-stop line] [-skip nlines]" << std::endl;
  std::cerr << std::endl;
  exit(0);
}

int main (int argc, char **argv){


  int i;
  unsigned int j, k, m;
  std::vector<Histogram*> histos;
  std::vector<std::string> inps;
  std::string currArg;
  std::vector<std::string> s;
  std::vector<std::vector<double> > data;
  std::ifstream inpFile;
  std::istream *finp;
  std::vector<std::vector<unsigned int> > cols;
  std::vector<unsigned int> bins;
  std::vector<double> res;
  bool verbose;
  std::string delim;
  std::string line;
  bool separate;
  HistoFormatEnum format;
  double temp;
  unsigned int start;
  unsigned int stop;
  unsigned int skip;
  unsigned int nline;

  histos.clear();
  inps.clear();
  finp=NULL;
  cols.clear();
  verbose=false;
  delim=" \t";
  line.clear();
  separate=false;
  bins.clear();
  res.clear();
  format=COUNT;
  temp=300;
  start=0;
  stop=std::numeric_limits<unsigned int>::max();
  skip=0;
  nline=0;

  for (i=1; i<argc; i++){
    currArg=argv[i];
    if (currArg.compare("-h") == 0 || currArg.compare("-help") == 0){
      usage();
    }
    else if (currArg.compare("-bins") == 0){
      currArg=argv[++i];
      Misc::splitNum(currArg, ":", bins, false);
    }
    else if (currArg.compare("-res") == 0){
      //Resolution/Bin spacing
      Misc::splitNum(currArg, ":", res, false);
      bins.clear();
    }
    else if (currArg.compare("-delimiter") == 0 || currArg.compare("-delimit") == 0){
      currArg=argv[++i];
      delim=currArg;
    }
    else if (currArg.compare("-col") == 0 || currArg.compare("-cols") == 0){
      currArg=argv[++i];
      cols.resize(cols.size()+1); //Expand vector
      Misc::splitNum(currArg, ":", cols.at(cols.size()-1), false);
      for (j=0; j< cols.at(cols.size()-1).size(); j++){
        //Shift columns by one to match vector reference
        cols.at(cols.size()-1).at(j)--;
      }
      //Create histogram with one window and nDim = number of columns
      histos.push_back(new Histogram(1, cols.at(cols.size()-1).size()));
    }
    else if (currArg.compare("-verbose") == 0){
      currArg=argv[++i];
      verbose=true;
    }
    else if (currArg.compare("-separate") == 0){
      separate=true;
    }
    else if (currArg.compare("-prob") == 0 || currArg.compare("-probability") == 0){
      format=PROBABILITY;
    }
    else if (currArg.compare("-density") == 0){
      //Normalized probability density
      format=DENSITY;
    }
    else if (currArg.compare("-energy") == 0){
      format=ENERGY;
      currArg=argv[++i];
      if (Misc::isdouble(currArg)){
        std::stringstream(currArg) >> temp;
      }
      else{
        std::cerr << "Error: Unrecognized temperature \"";
        std::cerr << currArg << "\"!" << std::endl;
        usage();
      }
    }
    else if (currArg.compare("-start") == 0){
      currArg=argv[++i];
      std::stringstream(currArg) >> start;
    }
    else if (currArg.compare("-stop") == 0){
      currArg=argv[++i];
      std::stringstream(currArg) >> stop;
    }
    else if (currArg.compare("-skip") == 0){
      currArg=argv[++i];
      std::stringstream(currArg) >> skip;
    }
		else if (currArg.compare(0,1,"-") == 0){
      std::cerr << "Warning: Skipping unknown option \"" << currArg << "\"" << std::endl;
    }
    else{
      inps.push_back(currArg);
    }
  }

  data.resize(histos.size());
  for (k=0; k< histos.size(); k++){
    data.at(k).resize(cols.at(k).size());
  }

  //time_t begin;
  //double end;
  //time(&begin);
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
    nline=0;
    while (finp->good() && !(finp->eof())){
      getline(*finp, line);
      if (line.length() > 0){
        nline++;
        if (nline < start || nline > stop){
          continue;
        }
        else if (skip > 0 && nline % (skip+1) != 0){
          continue;
        }
        else{
          //Don't use splitNum, str->double conversion too costly
          Misc::splitStr(line, delim, s, false);
          for (k=0; k< histos.size(); k++){
            for (m=0; m< cols.at(k).size(); m++){
              if (cols.at(k).at(m) < s.size()){
                //Only convert necessary columns from string to double
                std::stringstream(s.at(cols.at(k).at(m))) >> data.at(k).at(m);
              }
              else{
                std::cerr << std::endl << "Error: Missing data in column ";
                std::cerr << cols.at(k).at(m) << std::endl;
                exit(1);
              }
            }
            //Append vector of data into to zeroth window of kth histogram
            histos.at(k)->appendData(data.at(k), 0);
          }
        }
      }
    }
    if (inps.at(j).compare("-") !=0 && inpFile.is_open()){
      inpFile.close();
    }

    if (separate == true){
      //Output histograms for this file and clear data
      for (k=0; k< histos.size(); k++){
        histos.at(k)->setBins(bins);
        histos.at(k)->genHisto();
        histos.at(k)->printHisto(format, temp);
        std::cout << std::endl << std::endl;
        delete histos.at(k);
        histos.at(k) = new Histogram (1, cols.at(k).size());
      }
    }
  }
  //end=difftime(time(0), begin);
  //std::cout << "Time = " << end << " seconds " << std::endl;

  if (separate == false){
    //Output histograms
    for (k=0; k< histos.size(); k++){
      histos.at(k)->setBins(bins);
      histos.at(k)->genHisto();
      histos.at(k)->printHisto(format, temp);
      std::cout << std::endl << std::endl;
    }
  }

  return 0;
}
