//Sean M. Law

#include "PDB.hpp"

#include <iostream>
#include <iomanip>

int PDB::writePDB (Molecule& mol, std::string *ofile){
  std::ostringstream out;
  Atom atm;

  out.clear();
  for (unsigned int i=0; i<mol.getAtmVecSize(); i++){
    atm=mol.getAtom(i);
    out << std::setw(6) << std::left << atm.getRecName() << std::endl; 
  }

  if (ofile != 0){

  }
  else{
    std::cout << out.str(); 
  }
  return 0;
}
