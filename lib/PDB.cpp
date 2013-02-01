//Sean M. Law

#include "PDB.hpp"

#include <iostream>
#include <iomanip>

std::string PDB::writePDBFormat (Molecule& mol){
  std::ostringstream out;
  Atom atm;
  std::string lastChain="+";

  out.clear();
  for (unsigned int i=0; i<mol.getAtmVecSize(); i++){
    atm=mol.getAtom(i);
    if (atm.getChainId() != lastChain){
      if (lastChain != "+"){
        out << "TER" << std::endl;
      }
      lastChain=atm.getChainId();
    }
    out << std::setw(6) << std::left << atm.getRecName(); 
    out << std::setw(5) << std::right << atm.getAtmNum();
    out << " ";
    out << std::setw(4) << std::left << atm.getAtmName();
    out << std::setw(1) << std::left << atm.getAlt();
    out << std::setw(3) << std::left << atm.getResName();
    out << " ";
    out << std::setw(1) << std::left << atm.getChainId();
    out << std::setw(4) << std::right << atm.getResId();
    out << std::setw(1) << std::left << atm.getICode();
    out << "   ";
    out << std::fixed; //For setting precision
    out << std::setw(8) << std::right << std::setprecision(3) << atm.getX();
    out << std::setw(8) << std::right << std::setprecision(3) << atm.getY();
    out << std::setw(8) << std::right << std::setprecision(3) << atm.getZ();
    out << std::setw(6) << std::right << std::setprecision(2) << atm.getOccu();
    out << std::setw(6) << std::right << std::setprecision(2) << atm.getBFac();
    out << "      ";
    out << std::setw(4) << std::left << atm.getSegId();
    out << std::endl;
  }
  if (mol.getAtmVecSize()){
    out << "TER" << std::endl << "END" << std::endl;
  }

  return out.str();
}

void PDB::readPDB(Molecule& mol, std::string ifile, int model){
  std::ifstream pdbFile;
  std::istream* inp;
  std::string line;
  int currModel=0;
  std::string lastChain="+";
  std::string lastRes; //Check string.size()
  Atom atmEntry;

  if (ifile == "-"){ //Input from pipe
    inp=&std::cin;
  }
  else{ //Input from file
    pdbFile.open((ifile).c_str());
    inp=&pdbFile;
  }

  while (!(inp->eof())){

    getline(*inp,line);

    if (line.size() > 6 && line.compare(0,6,"MODEL ")==0){
      std::stringstream(line.substr(10,4)) >> currModel;
      if (model==0){
        model=1; //Use first model if undefined
      }
    }
    else if (currModel==model && line.size() >= 54 && (line.compare(0,4,"ATOM")==0 || line.compare(0,6,"HETATM")==0)){
      atmEntry=PDB::processAtomLine(line);
      mol.addAtom(atmEntry);
      //mol.addChain(PDB::compareChains(lastChain,atmEntry.getChainId()));
      atmEntry.reset();
    }
    else{
      continue;
    }
  }

  if (ifile != "-"){
    pdbFile.close();
  }
}

Atom PDB::processAtomLine (std::string line){
  double x,y,z;
  std::string recname; //Record name: "ATOM  ", "HETATM"
  int  atmnum; //Atom serial number
  int  resid; //Residue sequence number
  double occu; //Occupancy
  double bfac; //B-factor or Temperature factor
  std::string segid; //Segment identifier
  Atom atmEntry;

  atmEntry.reset();

  //substr: first character is denoted by a value of 0 (not 1)
  atmEntry.setRecName(line.substr(0,6));
  std::stringstream(line.substr(6,5)) >> atmnum;
  atmEntry.setAtmNum(atmnum);
  atmEntry.setAtmName(line.substr(12,4));
  atmEntry.setAlt(line.substr(16,1));
  atmEntry.setResName(line.substr(17,3));
  atmEntry.setChainId(line.substr(21,1));
  std::stringstream(line.substr(22,4)) >> resid;
  atmEntry.setResId(resid);
  atmEntry.setICode(line.substr(26,1));
  std::stringstream(line.substr(30,8)) >> x;
  std::stringstream(line.substr(38,8)) >> y;
  std::stringstream(line.substr(46,8)) >> z;
  atmEntry.setCoor(Vector(x,y,z));
  if (line.size() >= 60){
    std::stringstream(line.substr(54,6)) >> occu;
    atmEntry.setOccu(occu);
  }
  if (line.size() >= 66){
    std::stringstream(line.substr(60,6)) >> bfac;
    atmEntry.setBFac(bfac);
  }
  if (line.size() >= 76){
    segid=line.substr(72,4);
    atmEntry.setSegId(segid);
  }
  atmEntry.setSel(1);

  return atmEntry;
}

/*
Chain PDB::compareChains(std::string ch1, std::string ch2){

}
*/
