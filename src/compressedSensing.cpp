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

#include "Molecule.hpp"
#include "Atom.hpp"
#include "Misc.hpp"
#include "Analyze.hpp"

#include <sstream>

void usage(){
  std::cerr << std::endl;
  std::cerr << "Usage:   compressedSensing [-options] <MOL2files>" << std::endl;
  std::cerr << "Options: [-bins min:max:incr]" << std::endl;
  std::cerr << "         [-prosel selection] [-ligsel selection]" << std::endl;
//  std::cerr << "         [-format type] [-chains]" << std::endl;
//  std::cerr << "         [-warnings]" << std::endl;
  std::cerr << std::endl;
  exit(0);
}

int main (int argc, char **argv){


  int i, b;
  unsigned int j,k;
  std::vector<std::string> mol2;
  std::string currArg;
  std::string prosel;
  std::string ligsel;
  std::string format;
//  bool chnFlag;
//  bool warnings;
  double min, max, width;
  int bins;
  std::vector<double> s;
  std::vector<std::string> atomTypes;
  std::map<std::string, unsigned int> histo; //Sparse histogram
  std::vector<std::vector<double> > pdist;
  std::stringstream angstrom;
  std::string key;

  mol2.clear();
  prosel=":protein.";
  ligsel=":^protein.";
//  chnFlag=false;
  format.clear();
  min=0.0;
  max=10.0;
  width=2.0;

  //There might be a cleaner way of doing this
  //but adding new atom types is easier/more obvious here
  atomTypes.push_back("C.3");
  atomTypes.push_back("C.2");
  atomTypes.push_back("C.1");
  atomTypes.push_back("C.ar");
  atomTypes.push_back("C.cat");
  atomTypes.push_back("N.3");
  atomTypes.push_back("N.2");
  atomTypes.push_back("N.1");
  atomTypes.push_back("N.ar");
  atomTypes.push_back("N.am");
  atomTypes.push_back("N.pl3");
  atomTypes.push_back("N.4");
  atomTypes.push_back("O.3");
  atomTypes.push_back("O.2");
  atomTypes.push_back("O.co2");
  atomTypes.push_back("S.3");
  atomTypes.push_back("S.2");
  atomTypes.push_back("S.O");
  atomTypes.push_back("S.O2");
  atomTypes.push_back("P.3");
  atomTypes.push_back("F");
  atomTypes.push_back("Cl");
  atomTypes.push_back("Br");
  atomTypes.push_back("I");

  for (i=1; i<argc; i++){
    currArg=argv[i];
    if (currArg.compare("-h") == 0 || currArg.compare("-help") == 0){
      usage();
    }
    else if (currArg.compare("-prosel") == 0){
      currArg=argv[++i];
      prosel=currArg;
    }
    else if (currArg.compare("-ligsel") == 0){
      currArg=argv[++i];
      ligsel=currArg;
    }
    else if (currArg.compare("-format") == 0){
      currArg=argv[++i];
      Misc::toupper(currArg);
      format=currArg;
    }
//    else if (currArg.compare("-chains") == 0){
//      chnFlag=true;
//    }
    else if (currArg.compare("-bins") == 0){
      currArg=argv[++i];
      Misc::splitNum(currArg, ":", s, false);
      if (s.size() == 3){
        if (s.at(0) < s.at(1)){
          min=s.at(0);
          max=s.at(1);
        }
        else{
          min=s.at(1);
          max=s.at(0);
        }
        width=s.at(2);
      }
      else{
        std::cerr << "Error: Skipping unknown \"-bins\" format\n" << std::endl;
      }
    }
    else if (currArg.compare(0,1,"-") == 0){
      std::cerr << "Warning: Skipping unknown option \"" << currArg << "\"" << std::endl;
    }
    else{
      mol2.push_back(currArg);
    }
  }

  if (mol2.size() == 0){
    std::cerr << std::endl << "Error: Please provide an input file" << std::endl << std::endl;
    usage();
  }

  bins=static_cast<int>((max-min)/width);

  Molecule *mol=new Molecule;

  for (j=0; j< mol2.size(); j++){
    mol->cat(Molecule::readMol2(mol2.at(j), format));
  }

  mol->assignAtmInx();
  Analyze::pairwiseDistance(mol, pdist);

  //mol->writePDB(chnFlag);

  mol->select(ligsel);
  Molecule *ligand=mol->copy(true);

  Molecule *protein=NULL;

  //Fill histogram
  for (b=bins; b > 0; b--){
    max=min+b*width; //max gets updated here!
    angstrom.str(""); //Clear stringstream
    angstrom << min+b*width;
    mol->select(ligsel+"~"+angstrom.str()+"&"+prosel);
    protein=mol->copy(true); //Contains protein within "max" angstroms of ligand
    for (j=0; j< ligand->getNAtom(); j++){
      for (k=0; k< protein->getNAtom(); k++){
        double dist=pdist.at(ligand->getAtom(j)->getAtmInx()).at(protein->getAtom(k)->getAtmInx());
        if (dist <= max && dist > max-width){
          key=ligand->getAtom(j)->getAtmType()+":";
          key=key+protein->getAtom(k)->getAtmType()+":";
          key=key+angstrom.str();
          if (histo.find(key) != histo.end()){
            histo.at(key)++;
          }
          else{
            histo.insert(std::pair<std::string, unsigned int>(key, 1));
          }
          /*
          std::cerr << max-width << " - " << max << " ";
          std::cerr << ligand->getAtom(j)->getSummary() << " ";
          std::cerr << protein->getAtom(k)->getSummary() << std::endl;
          */
        }
      }
    }
  }

  //Output feature vector
  for (b=bins; b > 0; b--){
    angstrom.str(""); //Clear stringstream
    angstrom << min+b*width;
    for (j=0; j< atomTypes.size(); j++){
      for (k=0; k< atomTypes.size(); k++){
        key=atomTypes.at(j)+":";
        key=key+atomTypes.at(k)+":";
        key=key+angstrom.str();
        if (histo.find(key) != histo.end()){
          std::cout << histo.at(key);
        }
        else{
          std::cout << "0";
        }
        if (b != 1 || j != atomTypes.size()-1 || k != atomTypes.size()-1){
          std::cout << ",";
        }
        else{
          std::cout << std::endl;
        }
      }
    }
  }

  delete protein;
  delete ligand;
  delete mol;

  return 0;
}
