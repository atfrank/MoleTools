//Sean M. Law

#include <Misc.hpp>
#include <Select.hpp>
#include <Molecule.hpp>

#include <sstream>
#include <algorithm>

void Select::makeSel (Molecule* mol, std::string selin){
  //Select *sel=new Select;
  unsigned int i, j;
  std::vector<std::string> nExpr;
  std::vector<std::string> all;
  std::vector<std::string> tmp;
  Selection expr; //Declare a selection expression
  int nExcl, nColon, nPeriod;
  bool negAll=false;

  //Split logical OR operator "_"
  nExpr=Misc::split(selin, "_");
  for (i=0; i< nExpr.size(); i++){

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
        negAll=true;
        nExpr.at(i)=nExpr.at(i).substr(1, std::string::npos);
      }
      else{
        std::cerr << "Warning: Ignoring negation (!) in middle of selection \"";
        std::cerr << nExpr.at(i) << std::endl;
      }
    }

    //Split Chains/Segments, Residues, and Atoms
    all=Misc::split(nExpr.at(i), ":.");

    //Chains
    std::vector<Atom *> cmpChain=Select::chainRDP(all.at(0), mol->getAtmVec());

    for (i=0; i< cmpChain.size(); i++){
      std::cerr << cmpChain.at(i)->getSummary() << std::endl;
    }
    //Residues
    //all.at(1)

    //Atoms
    //all.at(2)

    //Overall negation
    if (negAll == true){

    } 

  }

}

std::vector<Atom *> Select::chainRDP (const std::string &str, const std::vector<Atom *> &ref){
  std::vector<Atom *> cmpCurr, cmpNext;
  std::vector<Atom *> out(2*ref.size());
  std::vector<Atom *>::iterator it;
  std::string curr, next;
  size_t pos;
  unsigned int i;

  if (str.length() == 0){
    return ref;
  }
  else if ((pos=str.find("+")) != std::string::npos){
    curr=str.substr(0, pos);
    cmpCurr=Select::chainRDP(curr, ref);
    next=str.substr(pos+1, std::string::npos);
    cmpNext=Select::chainRDP(next, ref);
    it=std::set_union(cmpCurr.begin(), cmpCurr.end(), cmpNext.begin(), cmpNext.end(), out.begin());
    out.resize(it-out.begin());
    it=std::unique(out.begin(), out.end());
    out.resize(std::distance(out.begin(),it));
  }
  else if ((pos=str.find("/")) != std::string::npos){
    curr=str.substr(0, pos);
    cmpCurr=Select::chainRDP(curr, ref);
    next=str.substr(pos+1, std::string::npos);
    cmpNext=Select::chainRDP(next, ref);
    it=std::set_intersection(cmpCurr.begin(), cmpCurr.end(), cmpNext.begin(), cmpNext.end(), out.begin());
    out.resize(it-out.begin());
  }
  else if ((pos=str.find("^")) == 0){
    curr=str.substr(pos+1, std::string::npos);
    cmpCurr=Select::chainRDP(curr, ref);
    it=std::set_difference(ref.begin(), ref.end(), cmpCurr.begin(), cmpCurr.end(), out.begin());
    out.resize(it-out.begin());
  }
  else{
    out.clear();
    for (i=0; i< ref.size(); i++){
      if (str == ref.at(i)->getChainId()){
        out.push_back(ref.at(i));
      }
      else if (str == ref.at(i)->getSegId()){
        out.push_back(ref.at(i));
      }
      else{
        continue;
      }
    }
  }

  return out;
}

void Select::makeSelOld (Molecule* mol, std::string selin){
  Select *sel=new Select;
  unsigned int i, j, k;
  bool flag;

  sel->parseSelOld(selin); //Change this to parseSel if needed

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
        if (Misc::trim(atm->getChainId()) == p.chainids.at(k)){
          //Chain matches
          if (p.negChainIds.at(k) == false){
            //No negation
            flag=true;
            break;
          }
        }
        else{ 
          //Chain does not match
          if (p.negChainIds.at(k) == true){
            //Apply negation
            flag=true;
            break;
          }          
        }
      }
      if (flag == false && p.chainids.size()){
        if(p.negAll == true){
          atm->setSel(true);
        }
        continue; //Next atom
      }

      //Check segid
      flag=false;
      for (k=0; k< p.segids.size(); k++){
        if (Misc::trim(atm->getSegId()) == p.segids.at(k)){
          //Segment matches
          if (p.negSegIds.at(k) == false){
            //No negation
            flag=true;
            break;
          }
        }
        else{
          //Segment does not match
          if (p.negSegIds.at(k) == true){
            //Apply negation
            flag=true;
            break;
          }
        }
      }
      if (flag == false && p.segids.size()){
        if(p.negAll == true){
          atm->setSel(true);
        }
        continue; //Next atom
      }
     
      //Check resname
      flag=false;
      for (k=0; k< p.resnames.size(); k++){
        if (Misc::trim(atm->getResName()) == p.resnames.at(k)){
          if (p.negResNames.at(k) == false){
            flag=true;
            break;
          }
        }
        else{
          if (p.negResNames.at(k) == true){
            flag=true;
            break;
          }
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
          if (p.negResIds.at(k) == false){
            flag=true;
            break;
          }
        }
        else{
          if (p.negResIds.at(k) == true){
            flag=true;
            break;
          }
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
        if (Misc::trim(atm->getAtmName()) == p.atmnames.at(k)){
          if (p.negAtmNames.at(k) == false){
            flag=true;
            break;
          }
        }
        else{
          if (p.negAtmNames.at(k) == true){
            flag=true;
            break;
          }
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

void Select::parseSelOld (std::string selin){
  std::vector<std::string> nExpr;
  std::vector<std::string> all;
  std::vector<std::string> tmp;
  Selection expr; //Declare a selection expression
  unsigned int i, j;
  int nExcl, nColon, nPeriod;
  bool negtmp;

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
      negtmp=false;
      if (tmp.at(j).substr(0,1) == "^"){
        negtmp=true;
        tmp.at(j)=tmp.at(j).substr(1,std::string::npos);
      }

      if (tmp.at(j).length() == 0){
        continue;
      }
      else if (tmp.at(j).length() == 4){
        expr.negSegIds.push_back(negtmp);
        expr.segids.push_back(tmp.at(j));
      }
      else if (tmp.at(j).length() == 1){
        expr.negChainIds.push_back(negtmp);
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
      negtmp=false;
      if (tmp.at(j).substr(0,1) == "^"){
        negtmp=true;
        tmp.at(j)=tmp.at(j).substr(1,std::string::npos);
      }
      
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
            expr.negResIds.push_back(negtmp);
            expr.resids.push_back(resid);
          }
        }
      }
      else if (Misc::isdigit(tmp.at(j))){
        std::stringstream(tmp.at(j)) >> resid;
        expr.negResIds.push_back(negtmp);
        expr.resids.push_back(resid);
      } 
      else if (tmp.at(j).length() == 3){  //4 character resnames??
        expr.negResIds.push_back(negtmp);
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
      negtmp=false;
      if (tmp.at(j).substr(0,1) == "^"){
        negtmp=true;
        tmp.at(j)=tmp.at(j).substr(1,std::string::npos);
      }
      if (tmp.at(j).length() ==0){
        continue;
      }
      else if (tmp.at(j).length() <= 4 && tmp.at(j).length() > 0){
        expr.negAtmNames.push_back(negtmp);
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
