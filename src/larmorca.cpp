//Sean M. Law
//Aaron T. Frank
    
/*
This file is part of MoleTools.

MoleTools is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MoleTools is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MoleTools.  If not, see <http://www.gnu.org/licenses/>.
Compiled Using: g++ larmorca.cpp -o larmorca -Wall -O3 -I$MoleTools/lib -lmoletools -L$MoleTools/lib -I$MoleTools */


#include "Molecule.hpp"
#include "Atom.hpp"
#include "Analyze.hpp"
#include "LARMORCA.hpp"
#include "Trajectory.hpp"

#include <iostream>
#include <ctime>
#include <fstream>

void usage(){
  std::cerr << "Usages: Skipping unknown option"  << std::endl;
  std::cerr << "Usage:   larmorca [-options] <PDBfile>" << std::endl;
  std::cerr << "Options: [-csfile CSfile]" << std::endl;
  std::cerr << "         [-trj TRAJfile]" << std::endl;
  std::cerr << "         [-skip frames] [-start frame] [-stop frame]" << std::endl;  
  std::cerr << "         [-identification ID]" << std::endl;
  std::cerr << std::endl;
  exit(0);
}

int main (int argc, char **argv){


  int i;
  unsigned int j;
  std::vector<std::string> pdbs;
  std::string currArg;
  std::vector<std::string> trajs;
  int start=0;
  int stop=std::numeric_limits<int>::max();
  int skip=0;
  bool startFlag=false;
  unsigned int itrj;
  std::ifstream trjin;
  Trajectory *ftrjin;
  unsigned int nframe;
  AnalyzeLarmorca *anin;
  bool verbose;
  std::string fchemshift;
  std::string identification;
  
  trajs.clear();
  ftrjin=NULL;
  nframe=0;
  anin=NULL;
  verbose=false;
  fchemshift="";
  identification="test";

  for (i=1; i<argc; i++){
    currArg=argv[i];
    if (currArg.compare("-h") == 0 || currArg.compare("-help") == 0){
      usage();
    }
    else if (currArg.compare("-trj") == 0 || currArg.compare("-traj") == 0){
      currArg=argv[++i];
      trajs.push_back(currArg);
    }
    else if (currArg.compare("-skip") == 0){
      currArg=argv[++i];
      std::stringstream(currArg) >> skip;
    }
    else if (currArg.compare("-start") == 0){
      currArg=argv[++i];
      std::stringstream(currArg) >> start;
      start--;
      startFlag=true;
    }
    else if (currArg.compare("-stop") == 0){
      currArg=argv[++i];
      std::stringstream(currArg) >> stop;
    }
    else if (currArg.compare("-csfile") == 0 || currArg.compare("-c") == 0 ){
      currArg=argv[++i];
      fchemshift=currArg;
    }  
    else if (currArg.compare("-identification") == 0 || currArg.compare("-i") == 0 ){
      currArg=argv[++i];
      identification=currArg;
    }
    else if (currArg.compare("-verbose") == 0){
      verbose=true;
    }
    else if (currArg.compare(0,1,"-") == 0){
      std::cerr << "Warning: Skipping unknown option \"" << currArg << "\"" << std::endl;
    }
    else{
      pdbs.push_back(currArg);
    }
  }

  if (pdbs.size() == 0){
    std::cerr << std::endl << "Error: Please provide an input file" << std::endl << std::endl;
    usage();
  }


  Molecule *mol=NULL;
  anin=new AnalyzeLarmorca;

  if (trajs.size() > 0){
    if (pdbs.size() > 1){
      std::cerr << std::endl << "Warning: Only the first PDB structure is used for trajectory analysis" << std::endl << std::endl;
    }
    /* Trajectory analysis */
    /* instantiate LARMORD */
    anin->addSel(":.CA");
    mol=Molecule::readPDB(pdbs.at(0));
    anin->preAnalysis(mol, "");
    
    /* Process trajectories */
    for (itrj=0; itrj< trajs.size(); itrj++){
      trjin.open(trajs.at(itrj).c_str(), std::ios::binary);
      if (trjin.is_open()){
        ftrjin=new Trajectory;
        ftrjin->setMolecule(mol);
        if (ftrjin->findFormat(trjin) == true){
          ftrjin->readHeader(trjin);
          if (skip > 0 && startFlag == false){
            start=skip;
          }
          /* Loop through desired frames */
          for (i=start; i< ftrjin->getNFrame() && i< stop; i=i+1+skip){
            if( ftrjin->readFrame(trjin, i) == false){
              std::cerr << "Warning: EOF found before the next frame could be read" << std::endl;
              break;
            }
            nframe++;
            anin->clearMol();
            anin->preAnalysis(mol, "");
            anin->runAnalysisTest(nframe,fchemshift,identification);
          }
        }
        else{
          std::cerr << "Warning: Skipping unknown trajectory format \"";
          std::cerr << trajs.at(itrj) << "\"" << std::endl;
        }
        if (ftrjin != NULL){
          delete ftrjin;
        }
      }
      trjin.close();
    }
  }
  else { 
    Molecule *mol=NULL;
      for (j=0; j< pdbs.size(); j++){
        mol=Molecule::readPDB(pdbs.at(j));
        anin->addSel(":.CA");
        anin->preAnalysis(mol, "");
        anin->runAnalysisTest(j+1,fchemshift,identification);
        delete mol;
      }

      if (anin != NULL){
        delete anin;
      }

      return 0;
    }
}