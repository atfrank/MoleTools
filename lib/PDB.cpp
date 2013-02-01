//Sean M. Law

#include "PDB.hpp"

#include <iostream>
#include <iomanip>

std::string PDB::writePDBFormat (Molecule& mol){
  std::ostringstream out;
  Atom atm;
  char lastChain='+';

  out.clear();
  for (unsigned int i=0; i<mol.getAtmVecSize(); i++){
    atm=mol.getAtom(i);
    if (atm.getChainId() != lastChain){
      if (lastChain != '+'){
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
