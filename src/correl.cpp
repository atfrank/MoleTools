//Sean M. Law

#include "Misc.hpp"

#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <cmath>
#include <limits>

void usage(){
  std::cerr << std::endl << std::endl;
  std::cerr << "Usage:   correl [-options] <file>" << std::endl;
  std::cerr << "Options: [-xcol column] [-ycol column]" << std::endl;
  std::cerr << "         [-start line] [-stop line]  [-skip lines]" << std::endl; 
	std::cerr << "         [-max value] [-min value]" << std::endl;
  exit(0);
}

int main (int argc, char **argv){


  int i;
  unsigned int j;
  unsigned int xcol;
  unsigned int ycol;
  unsigned int start, stop, skip;
	double max, min;
  std::string ifile;
  std::string currArg;
  std::ifstream inpFile;
  std::istream* inp;
  std::string line;
  std::vector<std::string> s;

  int n;
  unsigned int nline;
  std::vector<double> xvals;
  std::vector<double> yvals;
  double xval;
  double yval;
  double xavg;
  double yavg;
  double dx;
  double dy;
  double xx;
  double yy;
  double xy;
  double R;

  ifile.clear();
  xcol=0;
  ycol=1;
	start=0;
	stop=std::numeric_limits<unsigned int>::max();
  skip=0;
  n=0; //Counts the number datapoints used for correlation analysis
  nline=0; //Number of lines, used with skip
  xvals.clear();
  yvals.clear();
  xavg=0;
  yavg=0;
  xx=0;
  yy=0;
  xy=0;
	min=std::numeric_limits<double>::min();
	max=std::numeric_limits<double>::max();

  for (i=1; i<argc; i++){
    currArg=argv[i];
    if (currArg == "-h" || currArg == "-help"){
      usage();
    }
    else if (currArg == "-x" || currArg == "-xcol"){
      currArg=argv[++i];
      std::stringstream(currArg) >> xcol;
      xcol--;
    }
    else if (currArg == "-y" || currArg == "-ycol"){
      currArg=argv[++i];
      std::stringstream(currArg) >> ycol;
      ycol--;
    }
		else if (currArg == "-start"){
			currArg=argv[++i];
			std::stringstream(currArg) >> start;
		}
		else if (currArg == "-stop"){
			currArg=argv[++i];
			std::stringstream(currArg) >> stop;
		}
    else if (currArg == "-skip"){
      currArg=argv[++i];
      std::stringstream(currArg) >> skip;
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
		nline++;
    if (line.size() == 0 || nline < start || nline > stop){
      continue;
    }
		else if (skip > 0 && nline % (skip+1) != 0){
			continue;
		}
    else{
      Misc::splitStr(line, " \t", s, false); //Split on one or more consecutive whitespace
      if (s.size() <= xcol || s.size() <= ycol ){
        continue;
      }
      std::stringstream(s.at(xcol)) >> xval;
      std::stringstream(s.at(ycol)) >> yval;
			if (xval >= min && xval <= max && yval >= min && yval <= max){
      	xvals.push_back(xval);
      	yvals.push_back(yval);
      	xavg=xavg+xval;
      	yavg=yavg+yval;
      	n++;
			}
    }
  }

  if (ifile != "-"){
    inpFile.close();
  }

  xavg=xavg/n;
  yavg=yavg/n;
  
  for (j=0; j< xvals.size(); j++){
    dx=xvals.at(j)-xavg;
    dy=yvals.at(j)-yavg;
    xx=xx+dx*dx;
    yy=yy+dy*dy;
    xy=xy+dx*dy;
  }

  R=xy/(sqrt(xx*yy));

  std::cout << R << std::endl;
  
  return 0;
}
