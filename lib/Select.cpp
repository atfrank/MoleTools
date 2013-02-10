//Sean M. Law

#include <Misc.hpp>
#include <Select.hpp>
#include <Molecule.hpp>

void Select::makeSel (Molecule* mol, std::string sel){
  parseSel(sel);
}

void Select::parseSel (std::string sel){
  /*
  std::vector<std::string> qOr; //Logical OR
  std::vector<std::string> chnSeg;
  std::vector<std::string> chainId;
  std::vector<std::string> segId;
  std::vector<std::string> resName;
  std::vector<int> resId;
  std::vector<std::string> atmName;
  unsigned int i,j,k;

  qOr=Misc::split(sel, "_");
  for (i=0; i< qOr.size(); i++){
//    std::cout<< qOr.at(i) << std::endl;
    chnSeg=Misc::split(qOr.at(i), ":");
    for (j=0; j< chnSeg.size(); j++){
      
    }
  }
*/
  
}
