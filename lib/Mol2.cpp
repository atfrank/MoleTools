//Sean M. Law

#include "Mol2.hpp"

Mol2::Mol2(){
  chnMap.clear();
}

void Mol2::writeMol2Format (Molecule* mol, std::ostringstream &out, bool selFlag){
  Chain *chn;
  Residue *res;
  Atom *atm;
  std::string lastChain="+";
  int natom=0;
  int catom=0;

  out.clear();

  for (unsigned int i=0; i< mol->getChnVecSize(); i++){
    chn=mol->getChain(i);
    catom=0;
    //if(!chn->getSel()){continue;}
    for (unsigned int j=0; j< chn->getResVecSize(); j++){
      res=chn->getResidue(j);
      //if(!res->getSel()){continue;}
      for (unsigned int k=0; k< res->getAtmVecSize(); k++){
        atm=res->getAtom(k);
        if (selFlag == true && atm->getSel() == false){
					continue;
				}
				if (atm->getAtmNum() <= 99999){
					out << std::setw(6) << std::left << atm->getRecName();
        	out << std::setw(5) << std::right << atm->getAtmNum(); //Mol2 format
        	out << " "; //Mol2 format
				}
				else{
					out << std::setw(5) << std::left << atm->getRecName();
					out << std::setw(6) << std::right << atm->getAtmNum();
					out << " ";
				}
        out << std::setw(4) << std::left << atm->getAtmName();
        out << std::setw(1) << std::left << atm->getAlt();
				if (atm->getResName().length() < 4){
        	out << std::setw(3) << std::left << atm->getResName(); //Mol2 format
        	out << " "; //Mol2 format
				}
				else{
					out << std::setw(4) << std::left << atm->getResName();
				}
        out << std::setw(1) << std::left << atm->getChainId();
				if (atm->getResId() <= 999){
        	out << std::setw(4) << std::right << atm->getResId(); //Mol2 format
        	out << std::setw(1) << std::left << atm->getICode(); //Mol2 format
        	out << "   "; //Mol2 format
				}
//				else if (atm->getResId() <= 9999){
//					out << std::setw(5) << std::right << atm->getResId();
//          out << "   ";
//				}
				else if (atm->getResId() <= 99999){
					out << std::setw(5) << std::right << atm->getResId();
					out << "   ";
				}
				else{
					out << std::setw(6) << std::left << atm->getResId();
				}
        out << std::fixed; //For setting precision
        out << std::setw(8) << std::right << std::setprecision(3) << atm->getX();
        out << std::setw(8) << std::right << std::setprecision(3) << atm->getY();
        out << std::setw(8) << std::right << std::setprecision(3) << atm->getZ();
        out << std::setw(6) << std::right << std::setprecision(2) << atm->getOccu();
        out << std::setw(6) << std::right << std::setprecision(2) << atm->getBFac();
        out << "      ";
        out << std::setw(4) << std::left << atm->getSegId();
        out << std::endl;
        natom++;
        catom++;
      }
    }
    if (catom>0){
      out << "TER" << std::endl;
    }
  }

  if (natom > 0){
    out << "END" << std::endl;
  }

}

Molecule* Mol2::readMol2(std::string ifile){
  std::ifstream mol2File;
  std::istream* inp;
  std::string line;
  Molecule *mol=new Molecule;
  Chain *chnEntry=new Chain;
  Residue *resEntry=new Residue;
  Atom *atmEntry; //Created in heap by processAtomLine
  Atom *lastAtom;
  Mol2 mol2;
	bool readAtoms;

  atmEntry=NULL;
  lastAtom=NULL;
	mol->setCopyFlag(false);
	readAtoms=false;

  if (ifile == "-"){ //Input from pipe
    inp=&std::cin;
  }
  else{ //Input from file
    mol2File.open((ifile).c_str());
    inp=&mol2File;
  }

  while (inp->good() && !(inp->eof())){

    getline(*inp,line);
		if (line.size() >= 14 && line.compare(0,13,"@<TRIPOS>ATOM")==0){
			readAtoms=true;
    }
		else if (line.compare(0,9,"@<TRIPOS>") == 0){
			readAtoms=false;
		}
    else if (readAtoms == true && line.size() >= 44 ){
      //Atom
      atmEntry=mol2.processAtomLine(line, lastAtom);
			if (atmEntry == NULL){
				continue;
			}
      mol->addAtom(atmEntry);

      //Residue/Chain
      if (lastAtom != NULL){
				if(atmEntry->getChainId() != lastAtom->getChainId()){
        	//Store last
        	chnEntry->addResidue(resEntry);
        	mol->addResidue(resEntry);
        	mol->addChain(chnEntry);
        	//Create new
        	chnEntry=new Chain;
        	resEntry=new Residue;
      	}
      	else if (lastAtom->getResId() != atmEntry->getResId()){
        	chnEntry->addResidue(resEntry);
        	mol->addResidue(resEntry);
        	resEntry=new Residue;
					if (lastAtom->getResName() == "NME" || lastAtom->getResName() == "AHE" || atmEntry->getResName() == "UNK"){
						mol->addChain(chnEntry);
						chnEntry=new Chain;
					}
				}
				else{
					//Do Nothing
				}
      }
      else{
				//Do nothing
      }

			resEntry->addAtom(atmEntry);
			chnEntry->addAtom(atmEntry);

      //Update for next atom
      lastAtom=atmEntry;
    }
    else{
      continue;
    }
  }
  //Add remaining residues and chains
  if (resEntry->getAtmVecSize() > 0){
    chnEntry->addResidue(resEntry);
    mol->addResidue(resEntry);
    mol->addChain(chnEntry);
  }
  
  if (ifile != "-"){
    mol2File.close();
  }

  if (mol->getAtmVecSize() == 0){
    std::cerr << std::endl << "Error: The Mol2 file \"" << ifile << "\" ";
    std::cerr << "did not contain any valid atoms" << std::endl << std::endl;
    exit(1);
  }

  return mol;
}

Atom* Mol2::processAtomLine (std::string line, Atom* lastAtom){
  double x,y,z;
  std::string recname; //Record name: "ATOM  ", "HETATM"
  int  atmnum; //Atom serial number
  std::string chainid;
  int  resid; //Residue sequence number
  std::string segid; //Segment identifier
  Atom *atmEntry=new Atom;
	std::vector<std::string> s;
	Select *sel=new Select;
	std::vector<Atom *> ref;
	Molecule *mol=new Molecule;
	std::vector<Atom *> atmSel;

	sel->initKeys(mol); //mol is not actually used, needed to initialize selection keys

	Misc::splitStr(Misc::trim(line), " \t", s, false);
	std::stringstream(s.at(0)) >> atmnum;
  atmEntry->setAtmNum(atmnum);
  atmEntry->setAtmName(s.at(1));
	std::stringstream(s.at(2)) >> x;
  std::stringstream(s.at(3)) >> y;
  std::stringstream(s.at(4)) >> z;
  atmEntry->setCoor(Vector(x,y,z));

	if (s.size() >= 6){
		std::stringstream(s.at(6)) >> resid;
		atmEntry->setResId(resid);
	}
	if (s.size() >= 7){
		if (s.at(7).size() > 4){
			atmEntry->setResName(s.at(7).substr(0,4));
		}
		ref.clear();
		ref.push_back(atmEntry);
		atmSel=sel->recursiveDescentParser(":PROTEIN+NUCLEIC+ION+METAL+SOLVENT+TERMINI.", ref);
		if (atmSel.size() == 0){
			atmEntry->setResName(s.at(7).substr(0,3));	
			ref.clear();
	    ref.push_back(atmEntry);
	    atmSel=sel->recursiveDescentParser(":PROTEIN+NUCLEIC+ION+METAL+SOLVENT+TERMINI.", ref);
			if (atmSel.size() == 0){
			  atmEntry->setResName("UNK");	
			}
		}
	}

  atmEntry->setChainId(" ");

	std::stringstream() >> resid;
  atmEntry->setResId(resid);

	atmEntry->setSegId("");
  atmEntry->setSel(true);

  std::stringstream ss;
  ss << resid;

  atmEntry->setSummary(chainid+":"+atmEntry->getResName()+ss.str()+"."+Misc::trim(atmEntry->getAtmName()));

  return atmEntry;
}
