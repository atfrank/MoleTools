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
  std::vector<std::string> all;
  std::vector<std::string> tmp;
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
    if (nColon != 1 || nPeriod != 1){
      std::cerr << std::endl;
      std::cerr << "Warning: Skipping invalid selection expression syntax in \"";
      std::cerr << nExpr.at(i) << "\"" << std::endl;
      std::cerr << "Valid selection expression syntax:" << std::endl;
      std::cerr << "[!][[^]segid|chainid[+|/...]]<:>[[^]resid|resname|range[+|/...]]<.>[[^]atmname|key[+|/...]]";
      std::cerr << std::endl << std::endl;
      continue; //Evaluate next expression
    }

    //Split Chains/Segments from Residues/Atoms
    all=Misc::split(nExpr.at(i), ":.");

    for (j=0; j< all.size(); j++){
//      std::cerr << all.at(j) << std::endl;
    }

    /*
    //Chains/Segments
    tmp=Misc::split(all.at(0), "+");

    for (j=0; j< tmp.size(); j++){
      if (tmp.at(j).length() >= 4){
        expr.segids.push_back(tmp.at(j));
      }
      else if (tmp.at(j).length() >= 1){
        expr.chainids.push_back(tmp.at(j));
      }
      else if (tmp.at(j) == ""){
        continue; //Colon in front
      }
      else{
        std::cerr << "Warning: Skipping unrecognized Chain/Segment format\"";
        std::cerr << tmp.at(j) << "\"" << std::endl;
      }
    }

    //Residues
    tmp=Misc::split(all.at(1), "+");

    for (j=0; j< tmp.size(); j++){
      if (tmp.at(j).length() >= 4){
//        expr.segids.push_back(tmp.at(j));
      }
      else if (tmp.at(j).length() >= 1){
//        expr.chainids.push_back(tmp.at(j));
      }
      else if (tmp.at(j) == ""){
        continue; //Colon in front
      }
      else{
        std::cerr << "Warning: Skipping unrecognized Residue format\"";
        std::cerr << tmp.at(j) << "\"" << std::endl;
      }
    }

/*    
    //Atoms
    tmp=Misc::split(all.at(2), "+");

    for (j=0; j< tmp.size(); j++){
      if (tmp.at(j).length() >= 4){
        expr.segids.push_back(tmp.at(j));
      }
      else if (tmp.at(j).length() >= 1){
        expr.chainids.push_back(tmp.at(j));
      }
      else if (tmp.at(j) == ""){
        continue; //Colon in front
      }
      else{
        std::cerr << "Warning: Skipping unrecognized Atom format\"";
        std::cerr << tmp.at(j) << "\"" << std::endl;
      }
    }
*/
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
