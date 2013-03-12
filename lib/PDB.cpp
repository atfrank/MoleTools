//Sean M. Law

#include "PDB.hpp"

PDB::PDB(){
  chnMap.clear();
}

std::string PDB::writePDBFormat (Molecule* mol, bool selFlag ){
  std::ostringstream out;
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
        	out << std::setw(5) << std::right << atm->getAtmNum(); //PDB format
        	out << " "; //PDB format
				}
				else{
					out << std::setw(5) << std::left << atm->getRecName();
					out << std::setw(6) << std::right << atm->getAtmNum();
					out << " ";
				}
        out << std::setw(4) << std::left << atm->getAtmName();
        out << std::setw(1) << std::left << atm->getAlt();
				if (atm->getResName().length() < 4){
        	out << std::setw(3) << std::left << atm->getResName(); //PDB format
        	out << " "; //PDB format
				}
				else{
					out << std::setw(4) << std::left << atm->getResName();
				}
        out << std::setw(1) << std::left << atm->getChainId();
				if (atm->getResId() <= 999){
        	out << std::setw(4) << std::right << atm->getResId(); //PDB format
        	out << std::setw(1) << std::left << atm->getICode(); //PDB format
        	out << "   "; //PDB format
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

  return out.str();
}

Molecule* PDB::readPDB(std::string ifile, int model){
  std::ifstream pdbFile;
  std::istream* inp;
  std::string line;
  int currModel=0;
  Molecule *mol=new Molecule;
  Chain *chnEntry=new Chain;
  Residue *resEntry=new Residue;
  Atom *atmEntry; //Created in heap by processAtomLine
  Atom *lastAtom;
  PDB pdb;

  atmEntry=NULL;
  lastAtom=NULL;
	mol->setCopyFlag(false);

  if (ifile == "-"){ //Input from pipe
    inp=&std::cin;
  }
  else{ //Input from file
    pdbFile.open((ifile).c_str());
    inp=&pdbFile;
  }

  while (inp->good() && !(inp->eof())){

    getline(*inp,line);
    if (line.size() > 6 && line.compare(0,6,"MODEL ")==0){
      std::stringstream(line.substr(10,4)) >> currModel;
      if (model==0){
        model=1; //Use first model if undefined
      }
    }
    else if (currModel==model && line.size() >= 54 && (line.compare(0,4,"ATOM")==0 || line.compare(0,6,"HETATM")==0 || line.compare(0,5,"HETAT")==0)){
      //Atom
      atmEntry=pdb.processAtomLine(line, lastAtom);
      mol->addAtom(atmEntry);

      //Residue/Chain
      if (lastAtom != NULL && atmEntry->getChainId() != lastAtom->getChainId()) {
        //Store last
        chnEntry->addResidue(resEntry);
        mol->addResidue(resEntry);
        mol->addChain(chnEntry);
        //Create new
        chnEntry=new Chain;
        resEntry=new Residue;
      }
      else if (lastAtom != NULL && lastAtom->getResId() != atmEntry->getResId()){
        chnEntry->addResidue(resEntry);
        mol->addResidue(resEntry);
        resEntry=new Residue;
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
    pdbFile.close();
  }

  if (mol->getAtmVecSize() == 0){
    std::cerr << std::endl << "Error: The PDB file \"" << ifile << "\" ";
    std::cerr << "did not contain any valid atoms" << std::endl << std::endl;
    exit(1);
  }

  return mol;
}

Atom* PDB::processAtomLine (std::string line, Atom* lastAtom){
  double x,y,z;
  std::string recname; //Record name: "ATOM  ", "HETATM"
  int  atmnum; //Atom serial number
  std::string chainid;
  int  resid; //Residue sequence number
  double occu; //Occupancy
  double bfac; //B-factor or Temperature factor
  std::string segid; //Segment identifier
  Atom *atmEntry=new Atom;

  //substr: first character is denoted by a value of 0 (not 1)
	if (Misc::isdigit(Misc::trim(line.substr(4,7)))){
		atmEntry->setRecName(line.substr(0,4));
		std::stringstream(line.substr(4,7)) >> atmnum;
	  atmEntry->setAtmNum(atmnum);
	}
	else{
  	atmEntry->setRecName(line.substr(0,6)); //PDB format
  	std::stringstream(line.substr(6,5)) >> atmnum; //PDB format
  	atmEntry->setAtmNum(atmnum); //PDB format
	}
  atmEntry->setAtmName(line.substr(12,4));
  atmEntry->setAlt(line.substr(16,1));
  //atmEntry->setResName(line.substr(17,3)); //PDB format
	atmEntry->setResName(Misc::trim(line.substr(17,4)));
  chainid=line.substr(21,1);
  if(lastAtom == NULL || lastAtom->getChainId() != chainid){
    //New chain, check if duplicate
    if (chnMap.find(chainid) == chnMap.end()){
      chnMap[chainid]=chnMap.size();
    }
    else{
      //Duplicate, likely HETATM ligand bound to same chain
      //Assign numeric chainid, this is valid PDB format
      //std::ostringstream ossChainId;
      //ossChainId<< chnMap.find(chainid)->second;
      //chainid=ossChainId.str();
      //Set chainid to blank
      chainid=" ";
    }
  }
  atmEntry->setChainId(chainid);
  atmEntry->setRealId(line.substr(21,1)); 
  //std::stringstream(line.substr(22,4)) >> resid; //PDB format
	std::stringstream(line.substr(22,6)) >> resid;
  atmEntry->setResId(resid);
  atmEntry->setICode(line.substr(26,1));
  std::stringstream(line.substr(30,8)) >> x;
  std::stringstream(line.substr(38,8)) >> y;
  std::stringstream(line.substr(46,8)) >> z;
  atmEntry->setCoor(Vector(x,y,z));
  if (line.size() >= 60){
    std::stringstream(line.substr(54,6)) >> occu;
    atmEntry->setOccu(occu);
  }
  if (line.size() >= 66){
    std::stringstream(line.substr(60,6)) >> bfac;
    atmEntry->setBFac(bfac);
  }
  if (line.size() >= 76){
    segid=line.substr(72,4);
    atmEntry->setSegId(segid);
  }
  atmEntry->setSel(true);
  std::stringstream ss;
  ss << resid;
  atmEntry->setSummary(chainid+":"+atmEntry->getResName()+ss.str()+"."+Misc::trim(atmEntry->getAtmName()));

  return atmEntry;
}
