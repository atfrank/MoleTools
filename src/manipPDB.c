//#include <cstdlib> //exit()
//#include <iostream>
#include <fstream>
#include <sstream> //string stream
#include <string>
using namespace std;

#include <pdb.h>

#define MAXINPFILES 4096
#define MAXLINESIZE 4096

void usage (){
  fprintf (stderr, "\nUsage:  manipPDB [options] <PDBfile>\n");
  fprintf (stderr, "Options: [-h || -help]\n");
  fprintf (stderr, "\n\n");
  return;
}

int main (int argc, char **argv){

  ifstream fpdb;
  int i;
  int model;
  string currArg;
  stringstream ss;
  
  for (i=1; i<argc; i++){
    currArg=argv[i];
    if (currArg == "-h" || currArg == "-help"){
      usage();
    }
    else if (currArg == "-model"){
      ss << argv[i++];
      ss >> model;
    }
    else{
      fpdb.open(currArg.c_str());
    }
  }
  
  if (!fpdb){
    fprintf (stderr, "HERE\n");
  }

  return 0;
}
