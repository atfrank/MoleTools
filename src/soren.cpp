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
#include "Residue.hpp"
#include "Atom.hpp"
#include "Misc.hpp"
#include "Analyze.hpp"
#include "Constants.hpp"

#include <sstream>
#include <fstream>
#include <ctime>

void usage(){
  std::cerr << std::endl;
  std::cerr << "Usage:   soren [-options] <PDBfile>" << std::endl;
  std::cerr << "Options: [-bins min:max:incr] [-norm]" << std::endl;
  std::cerr << "         [-sel selection]" << std::endl;
  std::cerr << "         [-add type1[:type2[...[typeN]]] | -use type1[:type2[...[typeN]]]" << std::endl;
  std::cerr << "         [-response value}" << std::endl;
  std::cerr << "         [-chains]" << std::endl;
  std::cerr << "         [-header]" << std::endl;
  std::cerr << "         [-pka COEFfile]" << std::endl;
  std::cerr << "         [-modelpka]" << std::endl;
  std::cerr << "         [-top file]" << std::endl;
  std::cerr << "         [-mol2]" << std::endl;
//  std::cerr << "         [-warnings]" << std::endl;
  std::cerr << std::endl;
  exit(0);
}

int main (int argc, char **argv){

  int i, b;
  unsigned int j,k;
  std::vector<std::string> pdbs;
  std::string currArg;
  std::string sel;
  std::string format;
  bool chnFlag;
//  bool warnings;
  double min, max, width;
  int bins;
  std::vector<double> s;
  std::vector<std::string> sadd;
  std::vector<std::string> atomTypes;
  std::map<std::string, unsigned int> histo; //Sparse histogram
  std::vector<std::vector<double> > pdist;
  std::stringstream angstrom;
  std::string key;
  bool responseFlag;
  double response;
  bool headerFlag;
  bool normFlag;
  double surfaceArea;
  std::string fcoef;
  std::ifstream coefFile;
  std::istream* coefinp;
  std::string line;
  std::vector<std::string> r;
  double tmp;
  std::map<std::string,double> coef;
  double prediction;
  bool modelFlag;
  std::map<std::string, std::string> model;
  std::string resName;
  double dist;
  std::string top;
  bool mol2Flag;

  pdbs.clear();
  sel=":.OD1+OD2+OE1+OE2+NZ+SG+ND1+ND2";
  format.clear();
  min=0.0;
  max=10.0;
  width=1.0;
  responseFlag=false;
  response=0.0;
  headerFlag=false;
  normFlag=false;
  fcoef.clear();
  coef.clear();
  modelFlag=false;
  model.clear();
  resName.clear();
  top.clear();
  mol2Flag=false;

  //There might be a cleaner way of doing this
  //but adding new atom types is easier/more obvious here.
  //Note that they are in alphanumeric order and sorted later!
  //atomTypes.push_back("Br");
  //atomTypes.push_back("C.1");
  atomTypes.push_back("C.2");
  atomTypes.push_back("C.3");
  atomTypes.push_back("C.ar");
  atomTypes.push_back("C.cat");
  //atomTypes.push_back("Cl");
  //atomTypes.push_back("F");
  //atomTypes.push_back("I");
  //atomTypes.push_back("N.1");
  //atomTypes.push_back("N.2");
  atomTypes.push_back("N.3");
  atomTypes.push_back("N.4");
  atomTypes.push_back("N.am");
  atomTypes.push_back("N.ar");
  atomTypes.push_back("N.pl3");
  atomTypes.push_back("O.2");
  atomTypes.push_back("O.3");
  atomTypes.push_back("O.co2");
  //atomTypes.push_back("P.3");
  //atomTypes.push_back("S.2");
  atomTypes.push_back("S.3");
  //atomTypes.push_back("S.O");
  //atomTypes.push_back("S.O2");

  model.insert(std::pair<std::string, std::string>("UNK","0,0,0,0,0"));
  model.insert(std::pair<std::string, std::string>("ASP","1,0,0,0,0"));
  model.insert(std::pair<std::string, std::string>("GLU","0,1,0,0,0"));
  model.insert(std::pair<std::string, std::string>("CYS","0,0,1,0,0"));
  model.insert(std::pair<std::string, std::string>("HIS","0,0,0,1,0"));
  model.insert(std::pair<std::string, std::string>("LYS","0,0,0,0,1"));

  for (i=1; i<argc; i++){
    currArg=argv[i];
    if (currArg.compare("-h") == 0 || currArg.compare("-help") == 0){
      usage();
    }
    else if (currArg.compare("-sel") == 0){
      currArg=argv[++i];
      sel=currArg;
    }
    else if (currArg.compare("-format") == 0){
      currArg=argv[++i];
      format=Misc::toupper(currArg);
    }
    else if (currArg.compare("-chains") == 0){
      chnFlag=true;
    }
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
    else if (currArg.compare("-add") == 0 || currArg.compare("-use") == 0){
      if (currArg.compare("-use") == 0){
        atomTypes.clear();
      }
      currArg=argv[++i];
      Misc::splitStr(currArg, ":", sadd, false);
      atomTypes.insert(atomTypes.end(), sadd.begin(), sadd.end());
    }
    else if (currArg.compare("-response") == 0){
      currArg=argv[++i];
      std::stringstream(currArg) >> response;
      responseFlag=true;
    }
    else if (currArg.compare("-header") == 0){
      headerFlag=true;
    }
    else if (currArg.compare("-norm") == 0){
      normFlag=true;
    }
    else if (currArg.compare("-pka") == 0 || currArg.compare("-pkas") == 0 || currArg.compare("-pKa") == 0 || currArg.compare("-pKas") == 0){
      currArg=argv[++i];
      fcoef=currArg;
    }
    else if (currArg.compare("-modelpka") == 0 || currArg.compare("-modelpKa") == 0){
      modelFlag=true;
    }
    else if (currArg.compare("-top") == 0){
      currArg=argv[++i];
      top=currArg;
    }
    else if (currArg.compare("-mol2") == 0){
      mol2Flag=true;
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

  bins=static_cast<int>((max-min)/width);

  //Read coefficients
  if (fcoef.length() > 0){
    headerFlag=false;
    responseFlag=false;
    coefFile.open(fcoef.c_str(), std::ios::in);
    coefinp=&coefFile;
    if (coefinp->good() == false){
      std::cerr << "Warning: there was a problem reading file \"" << fcoef << "\"" << std::endl;
      std::cerr << "Generating features instead" << std::endl;
      fcoef.clear();
    }
    while (coefinp->good() && !(coefinp->eof())){
      getline(*coefinp, line);
      Misc::splitStr(Misc::trim(line), ", \t", r, false);
      if (r.size() >= 2){
        std::stringstream(r.at(1)) >> tmp;
        coef.insert(std::pair<std::string, double>(r.at(0), tmp));
      }
    }
  }

  Molecule *mol=new Molecule;
 
/*
  std::clock_t start;
  double duration;
  start=std::clock();
*/

  for (j=0; j< pdbs.size(); j++){
    if (mol2Flag == true || pdbs.at(j).substr(pdbs.at(j).length()-4).compare("mol2") == 0){
      //Has ".mol2" extension
      //Only atoms are read, bonds are not read
      mol->cat(Molecule::readMol2(pdbs.at(j), format));
    }
    else{
      mol->cat(Molecule::readPDB(pdbs.at(j), format));
    }
  }

  if (chnFlag == true){
    mol->addMissingChainIds();
  }

  
  mol->select(":.^HYDROGEN");
  Molecule *cmol=mol->copy(true); //Copy without hydrogens

  if (top.length() > 0){
    cmol->readTopology(top, false);
    //atomTypes=cmol->getPrmtop().getAtomTypes(); //Get all atom types
    //atomTypes=cmol->getAtomTypes(); //Only get atom types in molecule
    atomTypes.clear();
    atomTypes.push_back("C");
    atomTypes.push_back("CA");
    atomTypes.push_back("CAI");
    atomTypes.push_back("CC");
    atomTypes.push_back("CP1");
    atomTypes.push_back("CP2");
		atomTypes.push_back("CP3");
		atomTypes.push_back("CPH1");
		atomTypes.push_back("CPH2");
		atomTypes.push_back("CPT");
		atomTypes.push_back("CT1");
		atomTypes.push_back("CT2");
		atomTypes.push_back("CT2A");
		atomTypes.push_back("CT3");
		atomTypes.push_back("CY");
		atomTypes.push_back("N");
		atomTypes.push_back("NC2");
		atomTypes.push_back("NH1");
		atomTypes.push_back("NH2");
		atomTypes.push_back("NH3");
		atomTypes.push_back("NR2");
		atomTypes.push_back("NY");
		atomTypes.push_back("O");
		atomTypes.push_back("OC");
		atomTypes.push_back("OH1");
		atomTypes.push_back("S");
  }

  //Sort, remove duplicates, and resize atom type vector
  std::sort(atomTypes.begin(), atomTypes.end());
  std::vector<std::string>::iterator it=std::unique(atomTypes.begin(), atomTypes.end());
  atomTypes.resize(std::distance(atomTypes.begin(),it));

  cmol->select(sel);
  Molecule *tmol=cmol->copy(true); //Titratable

  Residue *res=tmol->getResidue(0); //Only look at first residue
  resName=Misc::toupper(res->getResName());

  //Fill histogram
  for (j=0; j< res->getAtmVecSize(); j++){
    for (k=0; k< cmol->getAtmVecSize(); k++){
      if (std::find(res->getAtmVec().begin(), res->getAtmVec().end(), cmol->getAtom(k)) != res->getAtmVec().end()){
        //Skip titratable atoms
        continue;
      }
      dist=Analyze::distance(res->getAtom(j)->getCoor(), cmol->getAtom(k)->getCoor());
      if (dist >= min && dist <= max){
        if (dist != max){
          b=static_cast<int>((dist-min)/width)+1; //Bin right edge
        }
        else{
          b=bins; //Last bin
        }
        angstrom.str(""); //Clear stringstream
        angstrom << min+(b-1)*width; //Bin left edge
        angstrom << "-";
        angstrom << min+b*width;
        key=res->getAtom(j)->getAtmType()+":";
        key=key+cmol->getAtom(k)->getAtmType()+":";
        key=key+angstrom.str();
        if (histo.find(key) != histo.end()){
          histo.at(key)++;
        }
        else{
          histo.insert(std::pair<std::string, unsigned int>(key, 1));
        }
      }
    }
  }

/*
  duration = (std::clock() - start) / (double) CLOCKS_PER_SEC;
  std::cerr << "* " << cmol->getNAtom() << " " << duration << std::endl;
*/

  //Output header vector
  if (headerFlag == true){
    if (responseFlag == true){
      std::cout << "pKa,";
    }
    for (b=bins; b > 0; b--){
      angstrom.str(""); //Clear stringstream
      angstrom << min+(b-1)*width; //Bin left edge
      angstrom << "-";
      angstrom << min+b*width; //Bin right edge
      for (j=0; j< atomTypes.size(); j++){
        for (k=0; k< atomTypes.size(); k++){
          key=atomTypes.at(j)+":";
          key=key+atomTypes.at(k)+":";
          key=key+angstrom.str();
          std::cout << key;
          if (b != 1 || j != atomTypes.size()-1 || k != atomTypes.size()-1){
            std::cout << ",";
          }
        }
      }  
    }
    if (modelFlag == true){
      std::cout << ",ASP,GLU,CYS,HIS,LYS";
    }
    std::cout << std::endl;
  }

  //Output feature vector
  prediction=0.0;
  if (responseFlag == true){
    std::cout << response << ",";
  }
  if (fcoef.length() > 0 && coef.find("pKa") != coef.end()){
    prediction+=coef.at("pKa");
  }
  for (b=bins; b > 0; b--){
    angstrom.str(""); //Clear stringstream
    angstrom << min+(b-1)*width; //Bin left edge
    angstrom << "-";
    angstrom << min+b*width;
    //Surface area for the middle of the bin (i.e., between bin left and bin right)
    surfaceArea=4*PI*(pow(min+(b-0.5)*width, 2.0));
    for (j=0; j< atomTypes.size(); j++){
      for (k=0; k< atomTypes.size(); k++){
        key=atomTypes.at(j)+":";
        key=key+atomTypes.at(k)+":";
        key=key+angstrom.str();
        if (fcoef.length() > 0){
          if (coef.find(key) != coef.end() && histo.find(key) != histo.end()){
            prediction+=histo.at(key)*coef.at(key);
          }
        }
        else{
          if (histo.find(key) != histo.end()){
            if (normFlag == true){
              //Normalize count by surface area
              std::cout << histo.at(key)/surfaceArea;
            }
            else{
              std::cout << histo.at(key);
            }
          }
          else{
            std::cout << "0";
          }
          if (b != 1 || j != atomTypes.size()-1 || k != atomTypes.size()-1){
            std::cout << ",";
          }
        }
      }
    }
  }
  if (modelFlag == true){
    if (model.find(resName) != model.end()){
      std::cout << "," << model.at(resName);
    }
    else{
      std::cout << "," << model.at("UNK");
    }
  }

  if (fcoef.length() > 0){
    std::cout << prediction;
  }  
  std::cout << std::endl;

  delete tmol;
  delete cmol;
  delete mol;

  return 0;
}
