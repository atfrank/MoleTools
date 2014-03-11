//Sean M. Law

#include "Select.hpp"

void Select::makeSel (Molecule* mol, std::string selin){

  std::vector<Atom *> ref;
  unsigned int i;

  //Convert selection to uppercase
  Misc::toupper(selin);

	//Initialize special selection keys
  Select *sel=new Select;
	sel->initKeys(mol);

  ref=mol->getAtmVecClone(); //Always make a clone of the pointers and sort it!
  std::sort(ref.begin(), ref.end());

  //Passing mol->getAtmVec() directly won't work
  //because it is not properly sorted!
  std::vector<Atom *> atmSel=sel->recursiveDescentParser(selin, ref);

	if (atmSel.size() == 0){
		std::cerr << std::endl << "Error: Selection \"";
		std::cerr << selin << "\" did not match any atoms";
		std::cerr << std::endl << std::endl;
		exit(1);
	}

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
    out.clear();
    next=str.substr(pos+1, std::string::npos);
    cmpNext=Select::recursiveDescentParser(next, ref, group);
    if (cmpNext.size() == 0){
			std::cerr << std::endl << "Warning: Selection \"";
	    std::cerr << next << "\" did not match any atoms";
	    std::cerr << std::endl << std::endl;
      return out;
    }
    std::sort(cmpNext.begin(), cmpNext.end());    

    curr=str.substr(0, pos);
    cmpCurr=Select::recursiveDescentParser(curr, ref, group);
    if (cmpCurr.size() == 0){
      return out;
    }
    std::sort(cmpCurr.begin(), cmpCurr.end());

    std::set_intersection(cmpCurr.begin(), cmpCurr.end(), cmpNext.begin(), cmpNext.end(), back_inserter(out));
  }
  else if ((pos=str.find("_")) != std::string::npos){
    //Logical OR between expressions: A:1-5.CA_B:10-15.CA
    next=str.substr(pos+1, std::string::npos);
    cmpNext=Select::recursiveDescentParser(next, ref, group);
		if (cmpNext.size() == 0){
      std::cerr << std::endl << "Warning: Selection \"";
      std::cerr << next << "\" did not match any atoms";
      std::cerr << std::endl << std::endl;
    }
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
    out.clear();
    next=str.substr(pos+1, std::string::npos);
    cmpNext=Select::recursiveDescentParser(next, ref, "atom");
    if (cmpNext.size() == 0){
      return out;
    }
    std::sort(cmpNext.begin(), cmpNext.end());

    if (pos > 0){
      curr=str.substr(0, pos);
      cmpCurr=Select::recursiveDescentParser(curr, ref, "residue");
      if (cmpCurr.size() == 0){
        return out;
      }
    }
    else{
      cmpCurr=ref;
    }
    std::sort(cmpCurr.begin(), cmpCurr.end());

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
    out.clear();
    next=str.substr(pos+1, std::string::npos);
    cmpNext=Select::recursiveDescentParser(next, ref, group);
    if (cmpNext.size() == 0){
      return out;
    }
    std::sort(cmpNext.begin(), cmpNext.end());

    curr=str.substr(0, pos);
    cmpCurr=Select::recursiveDescentParser(curr, ref, group);
    if (cmpCurr.size() == 0){
      return out;
    }
    std::sort(cmpCurr.begin(), cmpCurr.end());

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
	else if (selKeysRes.find(str) != selKeysRes.end() && group.compare("residue") == 0){
	  out.clear();
    if (selKeysRes[str].length() > 0){
		  out=Select::recursiveDescentParser(selKeysRes[str], ref, group);
    }
	}
	else if (selKeysAtm.find(str) != selKeysAtm.end() && group.compare("atom") == 0){
		out.clear();
    if (selKeysAtm[str].length() > 0){
		  out=Select::recursiveDescentParser(selKeysAtm[str], ref, group);
    }
	}
  else{
    out.clear();
    if (group.compare("chain") == 0){
			for (i=0; i< ref.size(); i++){
        if (str.compare(ref.at(i)->getChainId()) == 0){
          out.push_back(ref.at(i));
        }
        else if (str.compare(ref.at(i)->getSegId()) == 0){
          out.push_back(ref.at(i));
        }
        else{
          continue;
        }
			}
    }
    else if (group.compare("residue") == 0){
			for (i=0; i< ref.size(); i++){
        int resnum;
        std::stringstream(str) >> resnum;
        if (Misc::isdigit(str) && resnum == ref.at(i)->getResId()){
          out.push_back(ref.at(i));
        }
        else if (str.compare(ref.at(i)->getResName()) == 0){
          out.push_back(ref.at(i));
        }
				else if (str.compare(Misc::trim(ref.at(i)->getRecName())) == 0){
          out.push_back(ref.at(i));
        }
        else{
          continue;
        }
			}
    }
    else if (group.compare("atom") == 0){
			for (i=0; i< ref.size(); i++){
        if (str.compare(Misc::trim(ref.at(i)->getAtmName())) == 0){
          out.push_back(ref.at(i));
        }
        else{
          continue;
        }
			}
    }
    else{
      //Do Nothing
    }
  }

  return out;
}

std::string Select::getSelValue(const std::string &key){
	std::string foundKey;

	return foundKey;
}

void Select::initKeys(Molecule *mol){

	unsigned int i;
	std::string atmname;
	std::vector<std::string> heavy;
	std::vector<std::string> H;
	std::vector<std::string> O;
	std::vector<std::string> N;
	std::vector<std::string> C;
	std::vector<std::string> S;
	std::vector<std::string> P;

	//Protein/Peptide Backbone Atoms
	selKeysAtm["BACKBONE"]="C N O CA HA HA1 HA2 HN1 HN OT1 OT2 OXT HT1 HT2 HT3";

	//Nucleic Acid Backbone Atoms
	selKeysAtm["BACKBONE"]+=" P O1P O2P O5' O5* O3' O3* C3' C3* H3' H3* C4' C4* H4' H4* C5' C5* H5* H5' H5'' H3T H5T";

	//Protein Sidechain Atoms
	selKeysAtm["SIDECHAIN"]="CB CD CD1 CD2 CE CE1 CE2 CE3 CG CG1 CG2 CH2 CZ CZ2 CZ3 HB HB1 HB2 HB3 HD1 HD11 HD12 HD13 HD2 HD21 HD22 HD23 HD3 HE HE1 HE2 HE21 HE22 HE3 HG HG1 HG11 HG12 HG13 HG2 HG21 HG22 HG23 HH HH11 HH12 HH2 HH21 HH22 HZ HZ1 HZ2 HZ3 ND1 ND2 NE NE1 NE2 NH1 NH2 NZ OD1 OD2 OE1 OE2 OG OG1 OH SD SG";

	//Nucleic Sidechain Atoms
	selKeysAtm["SIDECHAIN"]+=" CB CD CD1 CD2 CE CE1 CE2 CE3 CG CG1 CG2 CH2 CZ CZ2 CZ3 HB HB1 HB2 HB3 HD1 HD11 HD12 HD13 HD2 HD21 HD22 HD23 HD3 HE HE1 HE2 HE21 HE22 HE3 HG HG1 HG11 HG12 HG13 HG2 HG21 HG22 HG23 HH HH11 HH12 HH2 HH21 HH22 HZ HZ1 HZ2 HZ3 ND1 ND2 NE NE1 NE2 NH1 NH2 NZ OD1 OD2 OE1 OE2 OG OG1 OH SD SG";

	//Sugar Atoms
	selKeysAtm["SUGAR"]="C1' C1* O4' O4* H1' H1* C2' C2* H2' H2'' H2* C3' C3* H3' H3* C4' C4* H4' H4* O3* O3'";

	//Base Atoms
	selKeysAtm["BASE"]="C2 C4 C5 C5M C6 C8 H1 H2 H21 H22 H3 H41 H42 H5 H51 H52 H53 H6 H61 H62 H8 N1 N2 N3 N4 N6 N7 N9 O2 O4 O6";

	selKeysAtm["HEAVY"]="";
	selKeysAtm["HYDROGEN"]="";
	selKeysAtm["OXYGEN"]="";
	selKeysAtm["NITROGEN"]="";
	selKeysAtm["CARBON"]="";
	selKeysAtm["SULFUR"]="";
	selKeysAtm["SULPHUR"]="";
	selKeysAtm["PHOSPHORUS"]="";
	selKeysAtm["PHOSPHOROUS"]="";
	for (i=0; i< mol->getAtmVecSize(); i++){
		atmname=Misc::trim(mol->getAtom(i)->getAtmName());
		if (Select::atom(atmname, "heavy", heavy)){
			heavy.push_back(atmname);
     	selKeysAtm["HEAVY"]+=atmname+" ";
		}
		else if (Select::atom(atmname, "H", H)){
			H.push_back(atmname);
     	selKeysAtm["HYDROGEN"]+=atmname+" ";
		}
		else if (Select::atom(atmname, "O", O)){
     	O.push_back(atmname);
     	selKeysAtm["OXYGEN"]+=atmname+" ";
   	}
		else if (Select::atom(atmname, "N", N)){
     	N.push_back(atmname);
     	selKeysAtm["NITROGEN"]+=atmname+" ";
   	}
		else if (Select::atom(atmname, "C", C)){
     	C.push_back(atmname);
     	selKeysAtm["CARBON"]+=atmname+" ";
   	}
		else if (Select::atom(atmname, "S", S)){
     	S.push_back(atmname);
			selKeysAtm["SULFUR"]+=atmname+" ";
     	selKeysAtm["SULPHUR"]+=atmname+" ";
   	}
		else if (Select::atom(atmname, "P", P)){
     	P.push_back(atmname);
     	selKeysAtm["PHOSPHORUS"]+=atmname+" ";
			selKeysAtm["PHOSPHOROUS"]+=atmname+" ";
   	}
		else{
			continue;
		}
	}
	
	selKeysAtm["HEAVY"]=Misc::trim(selKeysAtm["HEAVY"]);
	selKeysAtm["HYDROGEN"]=Misc::trim(selKeysAtm["HYDROGEN"]);
	selKeysAtm["OXYGEN"]=Misc::trim(selKeysAtm["OXYGEN"]);
  selKeysAtm["NITROGEN"]=Misc::trim(selKeysAtm["NITROGEN"]);
  selKeysAtm["CARBON"]=Misc::trim(selKeysAtm["CARBON"]);
  selKeysAtm["SULFUR"]=Misc::trim(selKeysAtm["SULFUR"]);
	selKeysAtm["SULPHUR"]=Misc::trim(selKeysAtm["SULPHUR"]);
  selKeysAtm["PHOSPHORUS"]=Misc::trim(selKeysAtm["PHOSPHORUS"]);
	selKeysAtm["PHOSPHOROUS"]=Misc::trim(selKeysAtm["PHOSPHOROUS"]);

  //std::cerr << selKeysAtm["HYDROGEN"] << std::endl;

  selKeysRes["HETERO"]="HETATM HETAT HETA";
	selKeysRes["PEPTIDE"]="ALA CYS CYX VAL LEU ILE ASP GLU GLY GLN ASN HSD HSE HSP HIE HID HIS PRO TRP MET SER THR PHE TYR LYS ARG";
	selKeysRes["PROTEIN"]=selKeysRes["PEPTIDE"];
	//At physiological pH
	selKeysRes["BASIC"]="ARG LYS";
	selKeysRes["ACIDIC"]="ASP GLU";
	selKeysRes["CHARGED"]="ARG LYS ASP GLU";
	selKeysRes["HYDROPHOBIC"]="ALA VAL LEU ILE PHE GLY PRO CYS CYX MET TRP";
	selKeysRes["POLAR"]="SER THR TYR ASN GLN";
	selKeysRes["NUCLEIC"]="ADE THY CYT GUA URA A T C G U DA DT DC DG DU";
	selKeysRes["PURINE"]="ADE GUA A G DA DG";
	selKeysRes["PYRIMIDINE"]="THY CYT URA T C U DT DC DU";
	selKeysRes["WATER"]="TIP3 TIP HOH SPC SPCE TIP4 TIP5";
	selKeysRes["METAL"]="ZN FE NI MN CU CO CA BE";
	selKeysRes["ION"]="MG NA CL K SOD CLA CLM NAP";
	selKeysRes["SOLVENT"]="MG NA CL K SOD CLA CLM NAP TIP3 TIP HOH SPC SPCE TIP4 TIP5";
	selKeysRes["TERMINI"]="ACE ACP AHE CT2 NME CT3 FOR";

	//Replace spaces with "+"
	std::replace(selKeysAtm["BACKBONE"].begin(), selKeysAtm["BACKBONE"].end(), ' ', '+');
	std::replace(selKeysAtm["SIDECHAIN"].begin(), selKeysAtm["SIDECHAIN"].end(), ' ', '+');
	std::replace(selKeysAtm["SUGAR"].begin(), selKeysAtm["SUGAR"].end(), ' ', '+');
 	std::replace(selKeysAtm["BASE"].begin(), selKeysAtm["BASE"].end(), ' ', '+');
	std::replace(selKeysAtm["HEAVY"].begin(), selKeysAtm["HEAVY"].end(), ' ', '+');
	std::replace(selKeysAtm["HYDROGEN"].begin(), selKeysAtm["HYDROGEN"].end(), ' ', '+');
	std::replace(selKeysAtm["OXYGEN"].begin(), selKeysAtm["OXYGEN"].end(), ' ', '+');
	std::replace(selKeysAtm["NITROGEN"].begin(), selKeysAtm["NITROGEN"].end(), ' ', '+');
	std::replace(selKeysAtm["CARBON"].begin(), selKeysAtm["CARBON"].end(), ' ', '+');
	std::replace(selKeysAtm["SULFUR"].begin(), selKeysAtm["SULFUR"].end(), ' ', '+');
	std::replace(selKeysAtm["SULPHUR"].begin(), selKeysAtm["SULPHUR"].end(), ' ', '+');
	std::replace(selKeysAtm["PHOSPHORUS"].begin(), selKeysAtm["PHOSPHORUS"].end(), ' ', '+');
	std::replace(selKeysAtm["PHOSPHOROUS"].begin(), selKeysAtm["PHOSPHOROUS"].end(), ' ', '+');

  std::replace(selKeysRes["HETERO"].begin(), selKeysRes["HETERO"].end(), ' ', '+');
	std::replace(selKeysRes["PEPTIDE"].begin(), selKeysRes["PEPTIDE"].end(), ' ', '+');
	std::replace(selKeysRes["PROTEIN"].begin(), selKeysRes["PROTEIN"].end(), ' ', '+');
	std::replace(selKeysRes["BASIC"].begin(), selKeysRes["BASIC"].end(), ' ', '+');
	std::replace(selKeysRes["ACIDIC"].begin(), selKeysRes["ACIDIC"].end(), ' ', '+');
	std::replace(selKeysRes["CHARGED"].begin(), selKeysRes["CHARGED"].end(), ' ', '+');
	std::replace(selKeysRes["HYDROPHOBIC"].begin(), selKeysRes["HYDROPHOBIC"].end(), ' ', '+');
	std::replace(selKeysRes["POLAR"].begin(), selKeysRes["POLAR"].end(), ' ', '+');
	std::replace(selKeysRes["NUCLEIC"].begin(), selKeysRes["NUCLEIC"].end(), ' ', '+');
	std::replace(selKeysRes["PURINE"].begin(), selKeysRes["PURINE"].end(), ' ', '+');
	std::replace(selKeysRes["PYRIMIDINE"].begin(), selKeysRes["PYRIMIDINE"].end(), ' ', '+');
	std::replace(selKeysRes["WATER"].begin(), selKeysRes["WATER"].end(), ' ', '+');
	std::replace(selKeysRes["METAL"].begin(), selKeysRes["METAL"].end(), ' ', '+');
	std::replace(selKeysRes["ION"].begin(), selKeysRes["ION"].end(), ' ', '+');
	std::replace(selKeysRes["SOLVENT"].begin(), selKeysRes["SOLVENT"].end(), ' ', '+');
	std::replace(selKeysRes["TERMINI"].begin(), selKeysRes["TERMINI"].end(), ' ', '+');
	
}

bool Select::atom(const std::string &str, std::string typein, const std::vector<std::string> &AtomVec){
	size_t pos;
	std::string atmType;

	pos=str.find_first_not_of("0123456789");
	atmType=str.substr(pos,1);

	if (typein.compare("heavy") == 0 && atmType.compare("H") != 0 && std::find(AtomVec.begin(),AtomVec.end(), str) == AtomVec.end()){
		return true;
	}
	else if (atmType.compare(typein) == 0 && std::find(AtomVec.begin(),AtomVec.end(), str) == AtomVec.end()){
		return true;
	}
	else{
		return false;
	}
}
