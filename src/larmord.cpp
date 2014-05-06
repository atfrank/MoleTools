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
*/

/* Compile Using: g++ larmord.cpp -o larmord -Wall -O3 -I$MoleTools/lib -lmoletools -L$MoleTools/lib -I$MoleTools */
#include "Molecule.hpp"
#include "Atom.hpp"
#include "Analyze.hpp"
#include "LARMORD.hpp"
#include "Trajectory.hpp"

#include <iostream>
#include <cstdlib>
#include <fstream>

void usage(){
  std::cerr << "Usages: Skipping unknown option"  << std::endl;
  std::cerr << "Usage:   larmord [-options] <PDBfile>" << std::endl;
  std::cerr << "Options: [-csfile CSfile]" << std::endl;
  std::cerr << "         [-parmfile PARAMfile]" << std::endl;
  std::cerr << "         [-trj TRAJfile]" << std::endl;
  std::cerr << "         [-skip frames] [-start frame] [-stop frame]" << std::endl;  
  std::cerr << std::endl;
  exit(0);
}

int main (int argc, char **argv){
  int i;
  unsigned int f;
  unsigned int ainx;
  std::stringstream resid;
  std::vector<std::string> pdbs;
  std::string currArg;
  std::string fchemshift;
  std::string fparmfile;
  std::string nucleus;
  std::string resname;
  std::string atomname;
  std::string key;
  std::string identification;

  std::vector<std::string> trajs;
  int start;
  int stop=std::numeric_limits<int>::max();
  int skip;
  bool startFlag=false;
  unsigned int itrj;
  std::ifstream trjin;
  Trajectory *ftrjin;
  unsigned int nframe;

  start=0;
  skip=0;
  nframe=0;
  
  double alpha;
  double dist;
  double cspred;
  double randcs;
  double expcs;

  std::vector<std::vector<double> > neighborDistances;

  Molecule *neighbormol;
  neighbormol=NULL;
  
  Atom *ai, *aj;  
  ai=NULL;
  aj=NULL;
  fchemshift="";
  fparmfile="";
  identification="None";
  
  LARMORD *larm;
  larm=NULL;

  pdbs.clear();

  for (i=1; i<argc; i++){
    currArg=argv[i];
    if (currArg.compare("-h") == 0 || currArg.compare("-help") == 0){
      usage();
    }
    else if (currArg.compare("-csfile") == 0 || currArg.compare("-c") == 0 ){
      currArg=argv[++i];
      fchemshift=currArg;
    }  
    else if (currArg.compare("-parmfile") == 0 || currArg.compare("-p") == 0 ){
      currArg=argv[++i];
      fparmfile=currArg;
    }
    else if (currArg.compare("-identification") == 0 || currArg.compare("-i") == 0 ){
      currArg=argv[++i];
      identification=currArg;
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
  if (trajs.size() > 0){
    if (pdbs.size() > 1){
      std::cerr << std::endl << "Warning: Only the first PDB structure is used for trajectory analysis" << std::endl << std::endl;
    }
    /* Trajectory analysis */
    /* instantiate LARMORD */
    mol=Molecule::readPDB(pdbs.at(0));
    mol->selAll();
    larm = new LARMORD(mol,fchemshift,fparmfile);
    
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
            /* get distance matrix */
            mol->assignAtmInx();
            mol->selAll();
            Analyze::pairwiseDistance(mol, neighborDistances);

            mol->select(":.HEAVY");
            neighbormol= mol->copy();
            for (unsigned int j=0; j< mol->getAtmVecSize(); j++){
              ai = mol->getAtom(j);
              nucleus = ai->getAtmName();
              if (larm->getShiftAtom(nucleus)==true){
                resname = ai->getResName();
                resid << ai->getResId();
                key = resid.str()+":"+nucleus;
                resid.str("");
                cspred = 0.0;
                expcs = larm->getExperimentalCS(key);
                //std::cerr << " I am here " << key << " " << expcs << std::endl;
                if(fchemshift.length() < 0 || expcs != 0.0){
                    key = resname+":"+nucleus;
                    randcs = larm->getRandomShift(key);
                    if(randcs > 0){
                      ainx = ai->getAtmInx();
                      for (unsigned int l=0; l < neighbormol->getAtmVecSize(); l++){
                        aj = neighbormol->getAtom(l);
                        if(ai!=aj){
                          resname = aj->getResName();
                          atomname = aj->getAtmName();
                          alpha = larm->getAlpha(nucleus+":"+aj->getResName()+":"+aj->getAtmName());
                          dist = neighborDistances.at(ainx).at(aj->getAtmInx());
                          cspred = cspred + alpha/dist/dist/dist;
                        }
                      }
                      cspred = cspred + randcs;
                      std::cout << 0 << " " << nframe << " " << ai->getResId() << " " << ai->getResName() << " " << nucleus << " " << cspred << " " << expcs << " " <<  randcs << " " << identification << std::endl;
                    }
                }
              }
            }
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
    /* instantiate LARMORD */
    for (f=0; f< pdbs.size(); f++){  
      mol=Molecule::readPDB(pdbs.at(f));
      larm = new LARMORD(mol,fchemshift,fparmfile);
      //std::cerr << "Processing file \"" << pdbs.at(f) << "..." << std::endl;
      /* get distance matrix */
      mol->assignAtmInx();
      mol->selAll();
      Analyze::pairwiseDistance(mol, neighborDistances);

      mol->select(":.HEAVY");
      neighbormol= mol->copy();
  
      for (unsigned int j=0; j< mol->getAtmVecSize(); j++){
        ai = mol->getAtom(j);
        nucleus = ai->getAtmName();
        if (larm->getShiftAtom(nucleus)==true){
          resname = ai->getResName();
          resid << ai->getResId();
          key = resid.str()+":"+nucleus;
          resid.str("");
          cspred = 0.0;
          expcs = larm->getExperimentalCS(key);
          //std::cerr << " I am here " << key << " " << expcs << std::endl;
          if(fchemshift.length() < 0 || expcs != 0.0){
              key = resname+":"+nucleus;
              randcs = larm->getRandomShift(key);
              if(randcs > 0){
                ainx = ai->getAtmInx();
                for (unsigned int l=0; l < neighbormol->getAtmVecSize(); l++){
                  aj = neighbormol->getAtom(l);
                  if(ai!=aj){
                    resname = aj->getResName();
                    atomname = aj->getAtmName();
                    alpha = larm->getAlpha(nucleus+":"+aj->getResName()+":"+aj->getAtmName());
                    dist = neighborDistances.at(ainx).at(aj->getAtmInx());
                    cspred = cspred + alpha/dist/dist/dist;
                  }
                }
                cspred = cspred + randcs;
                std::cout << 0 << " " << f+1 << " " << ai->getResId() << " " << ai->getResName() << " " << nucleus << " " << cspred << " " << expcs << " " <<  randcs << " " << identification << std::endl;
              }
          }
        }
      }
      delete mol;
    }
  }
  
  if(larm!=NULL)
    delete larm;
  return 0;
}
