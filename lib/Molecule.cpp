//Sean M. Law

#include "Molecule.hpp"
#include "PDB.hpp"

int Molecule::readPDB (std::string ifile, int model){

  std::ifstream pdbFile;
  std::string line;
  int currModel;
  double x,y,z;
  char lastChain;
  Atom *atmEntry;
  std::string recname; //Record name: "ATOM  ", "HETATM"
  int  atmnum; //Atom serial number
  int  resid; //Residue sequence number
  double occu; //Occupancy
  double bfac; //B-factor or Temperature factor
  std::string segid; //Segment identifier

  if (ifile == "-"){

  }
  else{
    pdbFile.open((ifile).c_str());
    if (pdbFile.is_open()){
      lastChain='+';
      while (pdbFile.good()){
        getline(pdbFile,line);
        //Check the model
        if (line.size() > 6 && line.compare(0,6,"MODEL ")==0){
          std::stringstream(line.substr(10,4)) >> currModel;
          if (model==0){
            model=1; //Use first model if undefined
          }
        }
        if (model && currModel != model){
          continue;
        }
        //Process atom entry
        if (line.size() >= 54 && (line.compare(0,4,"ATOM")==0 || line.compare(0,6,"HETATM")==0)){
          atmEntry=new Atom;
          //substr: first character is denoted by a value of 0 (not 1)
          if ((line.substr(21,1))[0] != lastChain){
            lastChain=line.substr(21,1)[0];
            //cout << "CHAIN " << lastChain << endl;
          }
          std::stringstream(line.substr(6,5)) >> atmnum;
          atmEntry->setAtmNum(atmnum);
          atmEntry->setAtmName(line.substr(12,4));
          atmEntry->setAlt((line.substr(16,1))[0]);
          atmEntry->setResName(line.substr(17,3));
          atmEntry->setChainId((line.substr(21,1))[0]);
          std::stringstream(line.substr(22,4)) >> resid;
          atmEntry->setResId(resid);
          atmEntry->setICode((line.substr(26,1))[0]);
          std::stringstream(line.substr(30,8)) >> x;
          std::stringstream(line.substr(38,8)) >> y;
          std::stringstream(line.substr(46,8)) >> z;
          atmEntry->setCoor(Vector(x,y,z));
          if (line.size() >= 60){
            std::stringstream(line.substr(54,6)) >> occu;
            atmEntry->setOccu(occu);
          }
          if (line.size() >= 66){
            std::stringstream(line.substr(60,6)) >> bfac;
            atmEntry->setBFac(bfac);
          }
          if (line.size() >= 76){
            segid=line.substr(72,4);
            atmEntry->setSegId(segid);
          }
          //cout << x << " " << y << " " << z << endl;
          atmVec.push_back(*atmEntry);
        }
      }
      /*
      for (unsigned int i=0; i<atmVec.size(); i++){
      cout << atmVec.at(i).getCoor().x() << "  ";
      cout << atmVec.at(i).getCoor().y() << "  ";
      cout << atmVec.at(i).getCoor().z() << endl;
      }
      */
      pdbFile.close();
    }
  }
  return 0;
}

int Molecule::writeMolecule(std::string format){

  std::string out;

  if (format.compare("pdb")==0){
    out=PDB::writePDBFormat(*this);
  }
  else{
    std::cerr << std::endl;
    std::cerr << "Error: Unrecognized output format \"" << format << "\"";
    std::cerr << std::endl;
    return 1;
  }

  std::cout << out;

  return 0;
}

Atom Molecule::getAtom(int atmnum){
  return atmVec.at(atmnum);
}

unsigned int Molecule::getAtmVecSize(){
  return atmVec.size();
}
