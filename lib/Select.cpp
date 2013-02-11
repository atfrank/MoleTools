//Sean M. Law

#include <Misc.hpp>
#include <Select.hpp>
#include <Molecule.hpp>

#include <sstream>

void Select::makeSel (Molecule* mol, std::string selin){
  Select *sel=new Select;
  unsigned int i, j, k;
  bool flag;

  sel->parseSel(selin);

  mol->deselAll();

  for (i=0; i< sel->selVec.size(); i++){
    Select::Selection p=sel->selVec.at(i);
    for (j=0; j< mol->getAtmVecSize(); j++){
      Atom *atm=mol->getAtom(j);
      if (atm->getSel() == true){
        continue; //Already satisfied previous selection
      }
     
      //Check chainid
      flag=false;
      for (k=0; k< p.chainids.size(); k++){
        if (atm->getChainId().find(p.chainids.at(k)) != std::string::npos){
          flag=true;
          break;
        }
      }
      if (flag == false && p.chainids.size()){
        continue; //Next atom
      }
     
      //Check resname
      flag=false;
      for (k=0; k< p.resnames.size(); k++){
        if (atm->getResName().find(p.resnames.at(k)) != std::string::npos){
          flag=true;
          break;
        }
      }
      if (flag == false && p.resnames.size()){
        continue; //Next atom
      }

      //Check resid
      flag=false;
      for (k=0; k< p.resids.size(); k++){
        if (atm->getResId() == p.resids.at(k)){
          flag=true;
          break;
        }
      }
      if (flag == false && p.resids.size()){
        continue; //Next atom
      }


      //Check atmname
      flag=false;
      for (k=0; k< p.atmnames.size(); k++){
        if (atm->getAtmName().find(p.atmnames.at(k)) != std::string::npos){
          flag=true;
          break;
        }
      }
      if (flag == false && p.atmnames.size()){
        continue; //Next atom
      }
      
      //Passed ALL checks
      atm->setSel(true);
    }
  }
}

void Select::parseSel (std::string selin){
  std::vector<std::string> nExpr;
  std::vector<std::string> all;
  std::vector<std::string> tmp;
  Selection expr; //Declare a selection expression
  unsigned int i, j;
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
      //std::cerr << "[!]"; //Negation
      std::cerr << "[[^]segid|chainid[+|/...]]<:>"; //Chains/Segments
      std::cerr << "[[^]resid|resname|range[+|/...]]<.>"; //Residues
      std::cerr << "[[^]atmname|key[+|/...]]"; //Atoms
      std::cerr << std::endl << std::endl;
      continue; //Evaluate next expression
    }

    //Split Chains/Segments, Residues, and Atoms
    all=Misc::split(nExpr.at(i), ":.");

    //Chains/Segments
    tmp=Misc::split(all.at(0), "+");

    for (j=0; j< tmp.size(); j++){
      if (tmp.at(j).length() == 0){
        continue;
      }
      else if (tmp.at(j).length() == 4 || tmp.at(j).length() == 5){
        expr.segids.push_back(tmp.at(j));
      }
      else if (tmp.at(j).length() == 1 || tmp.at(j).length() == 2){
        expr.chainids.push_back(tmp.at(j));
      }
      else{
        std::cerr << "Warning: Skipping unrecognized Chain/Segment format\"";
        std::cerr << tmp.at(j) << "\"" << std::endl;
      }
    }

    //Residues
    tmp=Misc::split(all.at(1), "+");
    int resid, start, stop;

    for (j=0; j< tmp.size(); j++){
      if (tmp.at(j).length() == 0){
        continue;
      }
      else if (tmp.at(j).length() == 3 || tmp.at(j).length() == 4){
        expr.resnames.push_back(tmp.at(j));
      }
      else if (Misc::isdigit(tmp.at(j))){
        std::stringstream(tmp.at(j)) >> resid;
        expr.resids.push_back(resid);
      }
      else if (Misc::isrange(tmp.at(j))){
        std::vector<std::string> range;
        range=Misc::split(tmp.at(j), "-");
        if (!Misc::isdigit(range.at(0)) || !Misc::isdigit(range.at(1)) || range.at(0).length() == 0 || range.at(1).length() == 0){
          std::cerr << "Warning: Skipping unrecognized Residue range\"";
          std::cerr << tmp.at(j) << "\"" << std::endl;
        }
        else{
          std::stringstream(range.at(0)) >> start;
          std::stringstream(range.at(1)) >> stop;
          if (start > stop){
            resid=start;
            start=stop;
            stop=resid;
          }
          for (resid=start; resid<=stop; resid++){
            expr.resids.push_back(resid);
          }
        }
      }
      else{
        std::cerr << "Warning: Skipping unrecognized Residue format\"";
        std::cerr << tmp.at(j) << "\"" << std::endl;
      }
    }

    
    //Atoms
    tmp=Misc::split(all.at(2), "+");

    for (j=0; j< tmp.size(); j++){
      if (tmp.at(j).length() ==0){
        continue;
      }
      else if (tmp.at(j).length() <= 4 && tmp.at(j).length() > 0){
        expr.atmnames.push_back(tmp.at(j));
      }
      else{
        std::cerr << "Warning: Skipping unrecognized Atom format\"";
        std::cerr << tmp.at(j) << "\"" << std::endl;
      }
    }

    selVec.push_back(expr);
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
