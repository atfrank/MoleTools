//Sean M. Law

#include <Misc.hpp>
#include <Select.hpp>
#include <Molecule.hpp>

#include <sstream>
#include <algorithm>
//#include <cctype>

void Select::makeSel (Molecule* mol, std::string selin){

  std::vector<Atom *> ref;
  unsigned int i;

/*
  for(i=0; i< str.length(); i++){
    if (Misc::isalpha(str.at(i))){
      str.at(i)=std::to
    }
  }
*/
//  std::transform(selin.begin(), selin.end(), selin.begin(),[](char c) {return (std::toupper(c); ); });

  ref=mol->getAtmVec(); //Always make a copy of the pointers and sort it!
  std::sort(ref.begin(), ref.end());

  //Passing mol->getAtmVec() directly won't work because it is not properly sorted!
  std::vector<Atom *> atmSel=Select::recursiveDescentParser(selin, ref);

  mol->deselAll();
  for (i=0; i< atmSel.size(); i++){
    atmSel.at(i)->setSel(true);
  }
}

std::vector<Atom *> Select::recursiveDescentParser (const std::string &str, const std::vector<Atom *> &ref, const std::string &group){
  std::vector<Atom *> cmpCurr, cmpNext;
  std::vector<Atom *> out;
  std::vector<Atom *>::iterator it;
  std::string curr, next, start, end;
  size_t pos;
  unsigned int i;

  //For memory efficiency, always parse "next" first!
  if (str.length() == 0){
    return ref;
  }
  else if ((pos=str.find("&")) != std::string::npos){
    //Logical AND between expressions: A:1-10.CA&A:5-15.CA = A:5-10.CA
    next=str.substr(pos+1, std::string::npos);
    cmpNext=Select::recursiveDescentParser(next, ref, group);
    std::sort(cmpNext.begin(), cmpNext.end());    

    curr=str.substr(0, pos);
    cmpCurr=Select::recursiveDescentParser(curr, ref, group);
    std::sort(cmpCurr.begin(), cmpCurr.end());

    out.clear();
    std::set_intersection(cmpCurr.begin(), cmpCurr.end(), cmpNext.begin(), cmpNext.end(), back_inserter(out));
  }
  else if ((pos=str.find("_")) != std::string::npos){
    //Logical OR between expressions: A:1-5.CA_B:10-15.CA
    next=str.substr(pos+1, std::string::npos);
    cmpNext=Select::recursiveDescentParser(next, ref, group);
    std::sort(cmpNext.begin(), cmpNext.end());

    curr=str.substr(0, pos);
    cmpCurr=Select::recursiveDescentParser(curr, ref, group);
    std::sort(cmpCurr.begin(), cmpCurr.end());

    out.resize(cmpCurr.size()+cmpNext.size());
    it=std::set_union(cmpCurr.begin(), cmpCurr.end(), cmpNext.begin(), cmpNext.end(), out.begin());
    out.resize(it-out.begin());

    it=std::unique(out.begin(), out.end());
    out.resize(std::distance(out.begin(),it));
  }
  else if ((pos=str.find("!")) == 0){
    //Expression negation: !A:1-5.CA
    curr=str.substr(pos+1, std::string::npos);
    cmpCurr=Select::recursiveDescentParser(curr, ref, group);
    std::sort(cmpCurr.begin(), cmpCurr.end());

    out.clear();
    std::set_difference(ref.begin(), ref.end(), cmpCurr.begin(), cmpCurr.end(), back_inserter(out));
  }
  else if ((pos=str.find(":")) != std::string::npos){
    //Logical AND between Chains and Residues: A:GLY.
    next=str.substr(pos+1, std::string::npos);
    cmpNext=Select::recursiveDescentParser(next, ref, group);
    std::sort(cmpNext.begin(), cmpNext.end());

    if (pos > 0){
      curr=str.substr(0, pos);
      cmpCurr=Select::recursiveDescentParser(curr, ref, "chain");
    }
    else{
      cmpCurr=ref;
    }
    std::sort(cmpCurr.begin(), cmpCurr.end());

    out.clear();
    std::set_intersection(cmpCurr.begin(), cmpCurr.end(), cmpNext.begin(), cmpNext.end(), back_inserter(out));
  }
  else if ((pos=str.find(".")) != std::string::npos){
    //Logical AND between Residues and Atoms: :GLY.CA
    next=str.substr(pos+1, std::string::npos);
    cmpNext=Select::recursiveDescentParser(next, ref, "atom");
    std::sort(cmpNext.begin(), cmpNext.end());

    if (pos > 0){
      curr=str.substr(0, pos);
      cmpCurr=Select::recursiveDescentParser(curr, ref, "residue");
    }
    else{
      cmpCurr=ref;
    }
    std::sort(cmpCurr.begin(), cmpCurr.end());

    out.clear();
    std::set_intersection(cmpCurr.begin(), cmpCurr.end(), cmpNext.begin(), cmpNext.end(), back_inserter(out));
  }
  else if ((pos=str.find("+")) != std::string::npos){
    //Logical OR between Chains, or between Residues, or between Atoms: 
    //A+B:., :GLY+ALA., :1+2+3., :.CA+CB+C+N
    next=str.substr(pos+1, std::string::npos);
    cmpNext=Select::recursiveDescentParser(next, ref, group);
    std::sort(cmpNext.begin(), cmpNext.end());

    curr=str.substr(0, pos);
    cmpCurr=Select::recursiveDescentParser(curr, ref, group);
    std::sort(cmpCurr.begin(), cmpCurr.end());

    out.resize(cmpCurr.size()+cmpNext.size());
    it=std::set_union(cmpCurr.begin(), cmpCurr.end(), cmpNext.begin(), cmpNext.end(), out.begin());
    out.resize(it-out.begin());

    it=std::unique(out.begin(), out.end());
    out.resize(std::distance(out.begin(),it));
  }
  else if ((pos=str.find("/")) != std::string::npos){
    //Logical AND between Chains, or between Residues, or between Atoms:
    //A/B:., :GLY/ALA., :1/2/3., :.CA/CB/C/N
    next=str.substr(pos+1, std::string::npos);
    cmpNext=Select::recursiveDescentParser(next, ref, group);
    std::sort(cmpNext.begin(), cmpNext.end());

    curr=str.substr(0, pos);
    cmpCurr=Select::recursiveDescentParser(curr, ref, group);
    std::sort(cmpCurr.begin(), cmpCurr.end());

    out.clear();
    std::set_intersection(cmpCurr.begin(), cmpCurr.end(), cmpNext.begin(), cmpNext.end(), back_inserter(out));
  }
  else if ((pos=str.find("-")) != std::string::npos){
    //Residue identifier range: :1-10.
    start=str.substr(0, pos);
    end=str.substr(pos+1, std::string::npos);
    if (Misc::isdigit(start) && Misc::isdigit(end)){
      out=Select::recursiveDescentParser(Misc::processRange(start, end), ref, group);
    }
    else{
      return ref;
    }
  }
  else if ((pos=str.find("^")) == 0){
    //Chain, Residue, or Atom negation: ^A:., :^GLY., :^2., :.^CA
    curr=str.substr(pos+1, std::string::npos);
    cmpCurr=Select::recursiveDescentParser(curr, ref, group);
    std::sort(cmpCurr.begin(), cmpCurr.end());

    out.clear();
    std::set_difference(ref.begin(), ref.end(), cmpCurr.begin(), cmpCurr.end(), back_inserter(out));
  }
  else{
    out.clear();
    for (i=0; i< ref.size(); i++){
      if (group == "chain"){
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
      else if (group == "residue"){
        int resnum;
        std::stringstream(str) >> resnum;
        if (Misc::isdigit(str) && resnum == ref.at(i)->getResId()){
          out.push_back(ref.at(i));
        }
        else if (str == ref.at(i)->getResName()){
          out.push_back(ref.at(i));
        }
        else{
          continue;
        }
      }
      else if (group == "atom"){
        if (str == Misc::trim(ref.at(i)->getAtmName())){
          out.push_back(ref.at(i));
        }
        else{
          continue;
        }
      }
      else{
        continue;
      }
    }
  }

  return out;
}
