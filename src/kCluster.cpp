//Sean M. Law

#include "Molecule.hpp"

void usage(){
  exit(0);
}

int main (int argc, char **argv){


  int i;
  std::string pdb;
  std::string currArg;
  std::string sel;

  pdb.clear();

  for (i=1; i<argc; i++){
    currArg=argv[i];
    if (currArg.compare("-h") == 0 || currArg.compare("-help") == 0){
      usage();
    }
    else if (currArg.compare("-nsel") == 0){
      currArg=argv[++i];
      sel=currArg;
    }
		else if (currArg.compare(0,1,"-") == 0){
      std::cerr << "Warning: Skipping unknown option \"" << currArg << "\"" << std::endl;
    }
    else{
      pdb=currArg;
    }
  }

  if (pdb.length() == 0){
    std::cerr << std::endl << "Error: Please provide an input file" << std::endl << std::endl;
    usage();
  }

  Molecule *mol=Molecule::readPDB(pdb);

  if (sel.length() >0){
    mol->select(sel);
  }

  Molecule *cmol=mol->clone();
  

  return 0;
}
