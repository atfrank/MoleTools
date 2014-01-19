//Sean M. Law

#include "Molecule.hpp"
#include "Analyze.hpp"

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cstdlib>
#include <ctime>

void usage(){
  std::cerr << "Usage:   pcasso [-options] <pdbFile>" << std::endl;
  std::cerr << "Options: [-predict | -features]" << std::endl;
//	std::cerr << "         [-dssp dsspFile]" << std::endl;
//	std::cerr << "         [-trial]" << std::endl;
  std::cerr << std::endl;
  exit(0);
}

int main (int argc, char **argv){


  int i;
	unsigned int j;
  std::vector<std::string> pdbs;
  std::string currArg;
  std::string dssp;
	PcassoOutEnum out;

  dssp.clear();
	out=PREDICT;

  for (i=1; i<argc; i++){
    currArg=argv[i];
    if (currArg.compare("-h") == 0 || currArg.compare("-help") == 0){
      usage();
    }
    else if (currArg.compare("-dssp") == 0){
      currArg=argv[++i];
      dssp=currArg;
    }
		else if (currArg.compare("-predict") == 0 || currArg.compare("-prediction") == 0){
			out=PREDICT;
		}
		else if (currArg.compare("-features") == 0 || currArg.compare("-feature") == 0){
			out=FEATURES;
		}
    else{
      pdbs.push_back(currArg);
    }
  }

  if (pdbs.size() == 0){
    std::cerr << std::endl << "Error: Please provide an input file" << std::endl << std::endl;
    usage();
  }

  if (dssp.length() > 0){
    if (pdbs.size() > 1){
      std::cerr << "Warning: \"-dssp\" option can only be used with a single PDB input" << std::endl;
      std::cerr << "Warning: Only the first PDB structure was processed" << std::endl;
    }
    if (pdbs.size() > 0){
      //Only process first structure!
      Molecule *mol=Molecule::readPDB(pdbs.at(0));
      std::cerr << "Processing file \"" << pdbs.at(0) << "..." << std::endl;
			mol->pcasso(dssp, out); //Makes temporary clone with C-alpha only, and analyzes it
      delete mol;
    }
  }
  else{
		//Placed here for efficiency; construct trees once instead of calling mol->pcasso()
		AnalyzePcasso* anin=new AnalyzePcasso;
		anin->addSel(":.CA");
		anin->setOutType(out);
		
	  for (j=0; j< pdbs.size(); j++){
  	  Molecule *mol=Molecule::readPDB(pdbs.at(j));
		  std::cerr << "Processing file \"" << pdbs.at(j) << "..." << std::endl;
				
			//std::clock_t start;
			//double duration;
			//start=std::clock();

			//mol->pcasso("", out); //Removed for efficiency; avoid re-constructing trees

			anin->clearMol();
			anin->preAnalysis(mol, "");

			anin->runAnalysis();

			//duration = (std::clock() - start) / (double) CLOCKS_PER_SEC;
			//std::cerr << duration << std::endl;

		  delete mol;
    }
		delete anin;
	}

  return 0;
}
