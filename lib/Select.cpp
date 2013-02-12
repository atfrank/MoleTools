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
        if(p.negAll == true){
          atm->setSel(true);
        }
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
        if(p.negAll == true){
          atm->setSel(true);
        }
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
        if(p.negAll == true){
          atm->setSel(true);
        }
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
        if(p.negAll == true){
          atm->setSel(true);
        }
        continue; //Next atom
      }
      
      //Passed ALL checks
      if(p.negAll == false){
          atm->setSel(true);
      }
      else{
        atm->setSel(false);
      }
    }
  }
}

void Select::parseSel (std::string selin){
  std::vector<std::string> nExpr;
  std::vector<std::string> all;
  std::vector<std::string> tmp;
  Selection expr; //Declare a selection expression
  unsigned int i, j;
  int nExcl, nColon, nPeriod;

  //Split logical OR operator "_"
  nExpr=Misc::split(selin, "_");
  for (i=0; i< nExpr.size(); i++){
    expr.clear();

    //Check syntax
    nExcl=0;
    nColon=0;
    nPeriod=0;
    for (j=0; j< nExpr.at(i).length(); j++){ 
      if (nExpr.at(i).substr(j,1) == "!"){
        nExcl++;
      }
      if (nExpr.at(i).substr(j,1) == ":"){
        nColon++;
      }
      if (nExpr.at(i).substr(j,1) == "."){
        nPeriod++;
      }
    }
    if (nExcl > 1 || nColon != 1 || nPeriod != 1){
      std::cerr << std::endl;
      std::cerr << "Warning: Ignoring invalid selection expression syntax ";
      if (nPeriod == 0){
        std::cerr << "(missing period) ";
      }
      if (nColon == 0){
        std::cerr << "(missing colon) ";
      }
      std::cerr << "in \"" << nExpr.at(i) << "\"" << std::endl;
      std::cerr << "Valid selection expression syntax:" << std::endl;
      std::cerr << "[!]"; //Negation
      std::cerr << "[[^]segid|chainid[+|/...]]<:>"; //Chains/Segments
      std::cerr << "[[^]resid|resname|range[+|/...]]<.>"; //Residues
      std::cerr << "[[^]atmname|key[+|/...]]"; //Atoms
      std::cerr << std::endl << std::endl;
      continue; //Evaluate next expression
    }

    //Handle negation (!)
    if (nExcl == 1){
      if(nExpr.at(i).substr(0,1) == "!"){
        expr.negAll=true;
        nExpr.at(i)=nExpr.at(i).substr(1, std::string::npos);
      }
      else{
        std::cerr << "Warning: Ignoring negation (!) in middle of selection \"";
        std::cerr << nExpr.at(i) << std::endl;
      }
    }

    //Split Chains/Segments, Residues, and Atoms
    all=Misc::split(nExpr.at(i), ":.");

    //Chains/Segments
    tmp=Misc::split(all.at(0), "+");

    for (j=0; j< tmp.size(); j++){
      if (tmp.at(j).length() == 0){
        continue;
      }
      else if (tmp.at(j).length() == 5 && tmp.at(j).substr(0,1) == "^"){
        expr.negSegIds.push_back(true);
        expr.segids.push_back(tmp.at(j).substr(1,std::string::npos));
      }
      else if (tmp.at(j).length() == 4){
        expr.negSegIds.push_back(false);
        expr.segids.push_back(tmp.at(j));
      }
      else if (tmp.at(j).length() == 2 && tmp.at(j).substr(0,1) == "^"){
        expr.negChainIds.push_back(true);
        expr.chainids.push_back(tmp.at(j).substr(1,std::string::npos));
      }
      else if (tmp.at(j).length() == 1){
        expr.negChainIds.push_back(false);
        expr.chainids.push_back(tmp.at(j));
      }
      else{
        std::cerr << "Warning: Ignoring unrecognized Chain/Segment format \"";
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
      else if (Misc::isrange(tmp.at(j))){
        std::vector<std::string> range;
        range=Misc::split(tmp.at(j), "-");
        if (!Misc::isdigit(range.at(0)) || !Misc::isdigit(range.at(1)) || range.at(0).length() == 0 || range.at(1).length() == 0){
          std::cerr << "Warning: Ignoring unrecognized Residue range \"";
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
      else if (Misc::isdigit(tmp.at(j))){
        std::stringstream(tmp.at(j)) >> resid;
        expr.resids.push_back(resid);
      } 
      else if (tmp.at(j).length() == 3){  //4 character resnames??
        expr.resnames.push_back(tmp.at(j));
      }
      else{
        std::cerr << "Warning: Ignoring unrecognized Residue format \"";
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
        std::cerr << "Warning: Skipping unrecognized Atom format \"";
        std::cerr << tmp.at(j) << "\"" << std::endl;
      }
    }

    selVec.push_back(expr);
  }
  
}

void Select::Selection::clear() {
  //Clear struct
  negAll=false;
  negChainIds.clear();
  chainids.clear();
  negSegIds.clear();
  segids.clear();
  negResNames.clear();
  resnames.clear();
  negResIds.clear();
  resids.clear();
  negAtmNames.clear();
  atmnames.clear();
}
