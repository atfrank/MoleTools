//Sean M. Law

#include <Misc.hpp>
#include <Select.hpp>
#include <Molecule.hpp>

void Select::makeSel (Molecule* mol, std::string selin){
  Select sel; 
  sel.parseSel(selin);
}

void Select::parseSel (std::string selin){
  std::vector<std::string> nExpr;
  std::vector<std::string> chnSeg_ResAtm;
  std::vector<std::string> chnSeg;
  std::vector<std::string> resAtm;
  std::vector<std::string> res;
  Selection expr; //Declare a selection expression
  unsigned int i, j, k;
  int nColon, nPeriod;

  //Split logical OR operator "_"
  nExpr=Misc::split(selin, "_");
  for (i=0; i< nExpr.size(); i++){
    expr.clear();

    //Check syntax
    nColon=0;
    nPeriod=0;
    for (j=0; j< nExpr.at(i).length(); j++){ 
      if (nExpr.at(i).substr(j,1) == ":"){
        nColon++;
      }
      if (nExpr.at(i).substr(j,1) == "."){
        nPeriod++;
      }
    }
    if (nColon > 1 || nPeriod > 1){
      std::cerr << "Warning: Skipping invalid selection syntax in \"";
      std::cerr << nExpr.at(i) << "\"" << std::endl;
      continue; //Evaluate next expression
    }

    //Split Chains/Segments from Residues/Atoms
    if (nColon){
      chnSeg_ResAtm=Misc::split(nExpr.at(i), ":");

      chnSeg=Misc::split(chnSeg_ResAtm.at(0), "+");

      for (j=0; j< chnSeg.size(); j++){
        if (chnSeg.at(j).length() == 4){
          expr.segids.push_back(chnSeg.at(j));
        }
        else if (chnSeg.at(j).length() == 1){
          expr.chainids.push_back(chnSeg.at(j));
        }
        else if (chnSeg.at(j) == ""){
          continue; //Colon in front
        }
        else{
          std::cerr << "Warning: Skipping invalid Chain/Segment \"";
          std::cerr << chnSeg.at(j) << "\"" << std::endl;
        }
      }

      if (chnSeg_ResAtm.size() == 2){

        resAtm=Misc::split(chnSeg_ResAtm.at(1), ".");

        res=Misc::split(resAtm.at(0), "+");
        for (k=0; k< resAtm.size(); k++){
          
        }
      }
    }
    else if (nPeriod){
      
      //std::cerr <<
    }
    //push to selVec;
  }
  
}

void Select::Selection::clear() {
  //Clear struct
  chainids.clear();
  segids.clear();
  resnames.clear();
  resids.clear();
  atmnames.clear();
}
