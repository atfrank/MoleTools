/* 

  Copyright University of Michigan.
  This file is part of the Larmor software suite and is made available under license.
  University of Michigan (UM) TECHtransfer: phone: 734-763-0614 email: techtransfer@umich.edu.
  
  Author: Aaron T. Frank
     
*/

#include "Molecule.hpp"
#include "Chain.hpp"
#include "Residue.hpp"
#include "Atom.hpp"
#include "Misc.hpp"
#include "Analyze.hpp"
#include "Trajectory.hpp"

#include <iostream>
#include <cstdlib>
#include <fstream>
#include <limits>
#include <sstream>
#include <math.h>
#include <algorithm>


void usage(){
  std::cerr << "====================================================" << std::endl;
  std::cerr << "CONTACTPDB  v1.00" << std::endl;
  std::cerr << "(c) 2014 Aaron T. Frank and University of Michigan." << std::endl;
  std::cerr << "====================================================" << std::endl;
  std::cerr << "Usage:   contactPDB [-options] <-sel (selection) -sel(selection) ...> <PDBfile>" << std::endl;
  std::cerr << "Options: [-cutoff distance]" << std::endl;
  std::cerr << "         [-buffer distance]" << std::endl;
  std::cerr << "         [-rij_mins file]" << std::endl;
  std::cerr << "         [-csfile file]" << std::endl;
  std::cerr << "         [-trj TRAJfile]" << std::endl;
  std::cerr << "         [-skip frames] [-start frame] [-stop frame]" << std::endl;  
  std::cerr << std::endl;
  exit(0);
}

std::string print_value(Chain *chn, int j){
	if (chn->getResidue(j) == NULL){
		return("TER");
	} else {
		return(chn->getResidue(j)->getResName());
	}
}

std::string print_key(Chain *chn, int j){
	std::stringstream resid;
	resid.str("");
	resid << chn->getResidue(j)->getResId();
	return(resid.str()+":"+chn->getResidue(j)->getResName()+":"+chn->getChainId());
}

void setup_nearset_neighbors(Molecule *molin, std::map<std::string,std::string> &nearset_neighbors){
  Chain *chn;
  std::string key, value;
	//loop over residue
  for ( int i=0; i< molin->getChnVecSize(); i++){
    chn=molin->getChain(i);
    for ( int j=0; j< chn->getResVecSize(); j++){
      key=print_key(chn,j);
      value=print_value(chn,j-1)+" "+print_value(chn,j+1);
      nearset_neighbors.insert(std::pair<std::string,std::string>(key,value));
      //std::cout << key << " " << value << std::endl;
    }
  }
}

void load_chem_shift_file(const std::string fcsfile, std::map<std::string,double> &cs_data){
	std::ifstream csFile;
	std::istream* csinp;
	std::string line, key;
	std::vector<std::string> s;	
	double exp_cs;
	if (fcsfile.length() > 0){
		csFile.open(fcsfile.c_str(), std::ios::in);
		csinp=&csFile;
		if (csFile.is_open()) {
			while (csinp->good() && !(csinp->eof())){
				getline(*csinp, line);
				Misc::splitStr(line, "     ", s, false);
				if (s.size() > 20){					
					//resid+resname+nucleus
					key = Misc::trim(s.at(5))+":"+Misc::trim(s.at(6))+":"+Misc::trim(s.at(7));
					exp_cs = atof(Misc::trim(s.at(10)).c_str());
					cs_data.insert(std::pair<std::string,double>(key,exp_cs));
					//std::cout << key << " " << exp_cs << std::endl;	
				}
			}
		}
		else {
			std::cerr << "Error: Unable to open file: " << fcsfile << std::endl;	
			exit(0);		                
		}
	}
}

void load_rij_mins_file(const std::string frij_mins, std::map<std::string,double> &rij_mins){
	std::ifstream rijFile;
	std::istream* rijinp;
	std::string line, key;
	std::vector<std::string> s;
	double rij_min;
	if (frij_mins.length() > 0){
		rijFile.open(frij_mins.c_str(), std::ios::in);
		rijinp=&rijFile;
		if (rijFile.is_open()) {
			while (rijinp->good() && !(rijinp->eof()))
			{
				getline(*rijinp, line);
				Misc::splitStr(line, " ", s, true);
				if (s.size() ==  3){
					key = Misc::trim(s.at(0))+":"+Misc::trim(s.at(1));				
					rij_min = atof(Misc::trim(s.at(2)).c_str());
					rij_mins.insert(std::pair<std::string,double>(key,rij_min));	
				}
			}
		}
		else {
			std::cerr << "Error: Unable to open file: " << frij_mins << std::endl;	
			exit(0);		                
		}
	}
}

void load_rij_mins(std::map<std::string,double> &rij_mins){
	rij_mins.insert(std::pair<std::string,double>("ALA:ALA",5.53016896037621));
	rij_mins.insert(std::pair<std::string,double>("ALA:ARG",7.99282937428693));
	rij_mins.insert(std::pair<std::string,double>("ALA:ASN",6.42209730753677));
	rij_mins.insert(std::pair<std::string,double>("ALA:ASP",6.24228731014858));
	rij_mins.insert(std::pair<std::string,double>("ALA:CYS",6.03474166562018));
	rij_mins.insert(std::pair<std::string,double>("ALA:GLN",6.97489093438761));
	rij_mins.insert(std::pair<std::string,double>("ALA:GLU",6.76216333226547));
	rij_mins.insert(std::pair<std::string,double>("ALA:GLY",5.0454));
	rij_mins.insert(std::pair<std::string,double>("ALA:HIS",6.90844179009027));
	rij_mins.insert(std::pair<std::string,double>("ALA:ILE",6.74398196945294));
	rij_mins.insert(std::pair<std::string,double>("ALA:LEU",6.96703072607167));
	rij_mins.insert(std::pair<std::string,double>("ALA:LYS",6.95798161799458));
	rij_mins.insert(std::pair<std::string,double>("ALA:MET",7.13581766175641));
	rij_mins.insert(std::pair<std::string,double>("ALA:PHE",7.59460226708356));
	rij_mins.insert(std::pair<std::string,double>("ALA:PRO",6.09947219927688));
	rij_mins.insert(std::pair<std::string,double>("ALA:SER",5.69802371965081));
	rij_mins.insert(std::pair<std::string,double>("ALA:THR",6.27260002241178));
	rij_mins.insert(std::pair<std::string,double>("ALA:TRP",7.88765024133371));
	rij_mins.insert(std::pair<std::string,double>("ALA:TYR",7.68141700198845));
	rij_mins.insert(std::pair<std::string,double>("ALA:VAL",6.35575403358502));
	rij_mins.insert(std::pair<std::string,double>("ARG:ARG",10.3096104426092));
	rij_mins.insert(std::pair<std::string,double>("ARG:ASN",8.80114994669552));
	rij_mins.insert(std::pair<std::string,double>("ARG:ASP",9.32150122201476));
	rij_mins.insert(std::pair<std::string,double>("ARG:CYS",8.02384900977305));
	rij_mins.insert(std::pair<std::string,double>("ARG:GLN",9.20925886720557));
	rij_mins.insert(std::pair<std::string,double>("ARG:GLU",9.87029143386778));
	rij_mins.insert(std::pair<std::string,double>("ARG:GLY",7.3309));
	rij_mins.insert(std::pair<std::string,double>("ARG:HIS",9.29208625825792));
	rij_mins.insert(std::pair<std::string,double>("ARG:ILE",8.83512238842503));
	rij_mins.insert(std::pair<std::string,double>("ARG:LEU",8.72937408187280));
	rij_mins.insert(std::pair<std::string,double>("ARG:LYS",9.87929757406129));
	rij_mins.insert(std::pair<std::string,double>("ARG:MET",8.95369667914543));
	rij_mins.insert(std::pair<std::string,double>("ARG:PHE",9.32081638426007));
	rij_mins.insert(std::pair<std::string,double>("ARG:PRO",8.40480453520069));
	rij_mins.insert(std::pair<std::string,double>("ARG:SER",8.41463875187521));
	rij_mins.insert(std::pair<std::string,double>("ARG:THR",8.56649260283177));
	rij_mins.insert(std::pair<std::string,double>("ARG:TRP",9.30192277103352));
	rij_mins.insert(std::pair<std::string,double>("ARG:TYR",9.74909719974754));
	rij_mins.insert(std::pair<std::string,double>("ARG:VAL",8.39777260160191));
	rij_mins.insert(std::pair<std::string,double>("ASN:ASN",7.30932590910709));
	rij_mins.insert(std::pair<std::string,double>("ASN:ASP",7.45754013791645));
	rij_mins.insert(std::pair<std::string,double>("ASN:CYS",6.65408403839377));
	rij_mins.insert(std::pair<std::string,double>("ASN:GLN",7.79065155001673));
	rij_mins.insert(std::pair<std::string,double>("ASN:GLU",7.95217053037599));
	rij_mins.insert(std::pair<std::string,double>("ASN:GLY",5.94145));
	rij_mins.insert(std::pair<std::string,double>("ASN:HIS",8.12488508154973));
	rij_mins.insert(std::pair<std::string,double>("ASN:ILE",7.33122875744647));
	rij_mins.insert(std::pair<std::string,double>("ASN:LEU",7.60844334816334));
	rij_mins.insert(std::pair<std::string,double>("ASN:LYS",8.37350005558603));
	rij_mins.insert(std::pair<std::string,double>("ASN:MET",8.0206248042581));
	rij_mins.insert(std::pair<std::string,double>("ASN:PHE",8.22827995888641));
	rij_mins.insert(std::pair<std::string,double>("ASN:PRO",6.82946201584029));
	rij_mins.insert(std::pair<std::string,double>("ASN:SER",6.82014585475383));
	rij_mins.insert(std::pair<std::string,double>("ASN:THR",7.00847184912274));
	rij_mins.insert(std::pair<std::string,double>("ASN:TRP",8.77722599638233));
	rij_mins.insert(std::pair<std::string,double>("ASN:TYR",8.70761417716907));
	rij_mins.insert(std::pair<std::string,double>("ASN:VAL",7.08839685661337));
	rij_mins.insert(std::pair<std::string,double>("ASP:ASP",7.52945237439578));
	rij_mins.insert(std::pair<std::string,double>("ASP:CYS",6.63503425874513));
	rij_mins.insert(std::pair<std::string,double>("ASP:GLN",8.06185713376395));
	rij_mins.insert(std::pair<std::string,double>("ASP:GLU",8.05800408673562));
	rij_mins.insert(std::pair<std::string,double>("ASP:GLY",5.962));
	rij_mins.insert(std::pair<std::string,double>("ASP:HIS",8.33517300429068));
	rij_mins.insert(std::pair<std::string,double>("ASP:ILE",7.31865614457308));
	rij_mins.insert(std::pair<std::string,double>("ASP:LEU",7.27005616075514));
	rij_mins.insert(std::pair<std::string,double>("ASP:LYS",8.76702667145956));
	rij_mins.insert(std::pair<std::string,double>("ASP:MET",7.65542025127410));
	rij_mins.insert(std::pair<std::string,double>("ASP:PHE",7.99549577620782));
	rij_mins.insert(std::pair<std::string,double>("ASP:PRO",6.72529167773547));
	rij_mins.insert(std::pair<std::string,double>("ASP:SER",6.67812665433557));
	rij_mins.insert(std::pair<std::string,double>("ASP:THR",6.80695633864552));
	rij_mins.insert(std::pair<std::string,double>("ASP:TRP",8.71583469760037));
	rij_mins.insert(std::pair<std::string,double>("ASP:TYR",8.98090931946629));
	rij_mins.insert(std::pair<std::string,double>("ASP:VAL",6.87270779794984));
	rij_mins.insert(std::pair<std::string,double>("CYS:CYS",5.92915313474635));
	rij_mins.insert(std::pair<std::string,double>("CYS:GLN",7.68855249687910));
	rij_mins.insert(std::pair<std::string,double>("CYS:GLU",7.09347238760586));
	rij_mins.insert(std::pair<std::string,double>("CYS:GLY",5.45815));
	rij_mins.insert(std::pair<std::string,double>("CYS:HIS",7.96839666757598));
	rij_mins.insert(std::pair<std::string,double>("CYS:ILE",7.27120492760452));
	rij_mins.insert(std::pair<std::string,double>("CYS:LEU",7.42502191233944));
	rij_mins.insert(std::pair<std::string,double>("CYS:LYS",7.57931680544716));
	rij_mins.insert(std::pair<std::string,double>("CYS:MET",7.58733103114290));
	rij_mins.insert(std::pair<std::string,double>("CYS:PHE",8.08425391765340));
	rij_mins.insert(std::pair<std::string,double>("CYS:PRO",6.64449117384391));
	rij_mins.insert(std::pair<std::string,double>("CYS:SER",6.37376670791375));
	rij_mins.insert(std::pair<std::string,double>("CYS:THR",6.54375634177523));
	rij_mins.insert(std::pair<std::string,double>("CYS:TRP",8.48020330636592));
	rij_mins.insert(std::pair<std::string,double>("CYS:TYR",7.99041153829002));
	rij_mins.insert(std::pair<std::string,double>("CYS:VAL",6.73630377587013));
	rij_mins.insert(std::pair<std::string,double>("GLN:GLN",8.73740880516343));
	rij_mins.insert(std::pair<std::string,double>("GLN:GLU",8.48664796863744));
	rij_mins.insert(std::pair<std::string,double>("GLN:GLY",6.5025));
	rij_mins.insert(std::pair<std::string,double>("GLN:HIS",8.65943938752173));
	rij_mins.insert(std::pair<std::string,double>("GLN:ILE",7.87792823291935));
	rij_mins.insert(std::pair<std::string,double>("GLN:LEU",8.12916633218305));
	rij_mins.insert(std::pair<std::string,double>("GLN:LYS",8.75162064637847));
	rij_mins.insert(std::pair<std::string,double>("GLN:MET",8.27745314935588));
	rij_mins.insert(std::pair<std::string,double>("GLN:PHE",8.65440969192088));
	rij_mins.insert(std::pair<std::string,double>("GLN:PRO",7.48845166126614));
	rij_mins.insert(std::pair<std::string,double>("GLN:SER",7.35405161984769));
	rij_mins.insert(std::pair<std::string,double>("GLN:THR",7.59704898130436));
	rij_mins.insert(std::pair<std::string,double>("GLN:TRP",8.91579716136401));
	rij_mins.insert(std::pair<std::string,double>("GLN:TYR",9.25254522163313));
	rij_mins.insert(std::pair<std::string,double>("GLN:VAL",7.59908903121547));
	rij_mins.insert(std::pair<std::string,double>("GLU:GLU",8.62613954494614));
	rij_mins.insert(std::pair<std::string,double>("GLU:GLY",6.47625));
	rij_mins.insert(std::pair<std::string,double>("GLU:HIS",8.90131657431987));
	rij_mins.insert(std::pair<std::string,double>("GLU:ILE",7.61182008266884));
	rij_mins.insert(std::pair<std::string,double>("GLU:LEU",7.86493129259806));
	rij_mins.insert(std::pair<std::string,double>("GLU:LYS",9.40041331384807));
	rij_mins.insert(std::pair<std::string,double>("GLU:MET",8.07416032861812));
	rij_mins.insert(std::pair<std::string,double>("GLU:PHE",8.38823401234979));
	rij_mins.insert(std::pair<std::string,double>("GLU:PRO",7.40097875309113));
	rij_mins.insert(std::pair<std::string,double>("GLU:SER",7.36278454835660));
	rij_mins.insert(std::pair<std::string,double>("GLU:THR",7.47345425232561));
	rij_mins.insert(std::pair<std::string,double>("GLU:TRP",9.10856157078892));
	rij_mins.insert(std::pair<std::string,double>("GLU:TYR",9.31009445344886));
	rij_mins.insert(std::pair<std::string,double>("GLU:VAL",7.35083468642770));
	rij_mins.insert(std::pair<std::string,double>("GLY:GLY",4.5417));
	rij_mins.insert(std::pair<std::string,double>("GLY:HIS",6.57985));
	rij_mins.insert(std::pair<std::string,double>("GLY:ILE",6.13975));
	rij_mins.insert(std::pair<std::string,double>("GLY:LEU",6.29895));
	rij_mins.insert(std::pair<std::string,double>("GLY:LYS",6.77155));
	rij_mins.insert(std::pair<std::string,double>("GLY:MET",6.4769));
	rij_mins.insert(std::pair<std::string,double>("GLY:PHE",6.8353));
	rij_mins.insert(std::pair<std::string,double>("GLY:PRO",5.5798));
	rij_mins.insert(std::pair<std::string,double>("GLY:SER",5.38));
	rij_mins.insert(std::pair<std::string,double>("GLY:THR",5.69115));
	rij_mins.insert(std::pair<std::string,double>("GLY:TRP",7.26675));
	rij_mins.insert(std::pair<std::string,double>("GLY:TYR",7.12165));
	rij_mins.insert(std::pair<std::string,double>("GLY:VAL",5.78365));
	rij_mins.insert(std::pair<std::string,double>("HIS:HIS",8.55717572308617));
	rij_mins.insert(std::pair<std::string,double>("HIS:ILE",8.00152500055914));
	rij_mins.insert(std::pair<std::string,double>("HIS:LEU",8.25738684363020));
	rij_mins.insert(std::pair<std::string,double>("HIS:LYS",8.73923819920487));
	rij_mins.insert(std::pair<std::string,double>("HIS:MET",8.38112223875202));
	rij_mins.insert(std::pair<std::string,double>("HIS:PHE",8.70493588933860));
	rij_mins.insert(std::pair<std::string,double>("HIS:PRO",7.31353208860045));
	rij_mins.insert(std::pair<std::string,double>("HIS:SER",7.49712673990864));
	rij_mins.insert(std::pair<std::string,double>("HIS:THR",7.69748543290461));
	rij_mins.insert(std::pair<std::string,double>("HIS:TRP",9.31693768131733));
	rij_mins.insert(std::pair<std::string,double>("HIS:TYR",9.15176528904577));
	rij_mins.insert(std::pair<std::string,double>("HIS:VAL",7.67077605403691));
	rij_mins.insert(std::pair<std::string,double>("ILE:ILE",7.93421991976067));
	rij_mins.insert(std::pair<std::string,double>("ILE:LEU",8.13629150574564));
	rij_mins.insert(std::pair<std::string,double>("ILE:LYS",7.90776366647095));
	rij_mins.insert(std::pair<std::string,double>("ILE:MET",8.30766159679873));
	rij_mins.insert(std::pair<std::string,double>("ILE:PHE",8.60517530507573));
	rij_mins.insert(std::pair<std::string,double>("ILE:PRO",7.26814354404592));
	rij_mins.insert(std::pair<std::string,double>("ILE:SER",6.82215102934463));
	rij_mins.insert(std::pair<std::string,double>("ILE:THR",7.40818376111115));
	rij_mins.insert(std::pair<std::string,double>("ILE:TRP",9.2349162472828));
	rij_mins.insert(std::pair<std::string,double>("ILE:TYR",8.58373576779798));
	rij_mins.insert(std::pair<std::string,double>("ILE:VAL",7.53153017093052));
	rij_mins.insert(std::pair<std::string,double>("LEU:LEU",8.23384298695426));
	rij_mins.insert(std::pair<std::string,double>("LEU:LYS",7.98743593838137));
	rij_mins.insert(std::pair<std::string,double>("LEU:MET",8.38193272486145));
	rij_mins.insert(std::pair<std::string,double>("LEU:PHE",8.77038561935594));
	rij_mins.insert(std::pair<std::string,double>("LEU:PRO",7.48830069717538));
	rij_mins.insert(std::pair<std::string,double>("LEU:SER",7.15891084289977));
	rij_mins.insert(std::pair<std::string,double>("LEU:THR",7.61412588144332));
	rij_mins.insert(std::pair<std::string,double>("LEU:TRP",9.23033933481777));
	rij_mins.insert(std::pair<std::string,double>("LEU:TYR",8.75772000243145));
	rij_mins.insert(std::pair<std::string,double>("LEU:VAL",7.77376253472725));
	rij_mins.insert(std::pair<std::string,double>("LYS:LYS",9.54178982590319));
	rij_mins.insert(std::pair<std::string,double>("LYS:MET",8.43967395811725));
	rij_mins.insert(std::pair<std::string,double>("LYS:PHE",8.46690628498818));
	rij_mins.insert(std::pair<std::string,double>("LYS:PRO",7.76680673916998));
	rij_mins.insert(std::pair<std::string,double>("LYS:SER",7.97159185739172));
	rij_mins.insert(std::pair<std::string,double>("LYS:THR",8.01631217486102));
	rij_mins.insert(std::pair<std::string,double>("LYS:TRP",8.83839808903178));
	rij_mins.insert(std::pair<std::string,double>("LYS:TYR",9.29901654905215));
	rij_mins.insert(std::pair<std::string,double>("LYS:VAL",7.53528936963561));
	rij_mins.insert(std::pair<std::string,double>("MET:MET",8.49031267437161));
	rij_mins.insert(std::pair<std::string,double>("MET:PHE",8.88887965347402));
	rij_mins.insert(std::pair<std::string,double>("MET:PRO",7.71351514468974));
	rij_mins.insert(std::pair<std::string,double>("MET:SER",7.12799532329204));
	rij_mins.insert(std::pair<std::string,double>("MET:THR",7.70429613364061));
	rij_mins.insert(std::pair<std::string,double>("MET:TRP",9.38719396473260));
	rij_mins.insert(std::pair<std::string,double>("MET:TYR",8.8439362576805));
	rij_mins.insert(std::pair<std::string,double>("MET:VAL",7.94308189480210));
	rij_mins.insert(std::pair<std::string,double>("PHE:PHE",9.21214765037181));
	rij_mins.insert(std::pair<std::string,double>("PHE:PRO",7.72313302939487));
	rij_mins.insert(std::pair<std::string,double>("PHE:SER",7.57892543204881));
	rij_mins.insert(std::pair<std::string,double>("PHE:THR",8.18940376849220));
	rij_mins.insert(std::pair<std::string,double>("PHE:TRP",9.97302825041026));
	rij_mins.insert(std::pair<std::string,double>("PHE:TYR",9.29039126520243));
	rij_mins.insert(std::pair<std::string,double>("PHE:VAL",8.44741364761614));
	rij_mins.insert(std::pair<std::string,double>("PRO:PRO",6.71260953770053));
	rij_mins.insert(std::pair<std::string,double>("PRO:SER",6.41988704894438));
	rij_mins.insert(std::pair<std::string,double>("PRO:THR",6.85390020808575));
	rij_mins.insert(std::pair<std::string,double>("PRO:TRP",8.09578489266612));
	rij_mins.insert(std::pair<std::string,double>("PRO:TYR",8.30378623251600));
	rij_mins.insert(std::pair<std::string,double>("PRO:VAL",6.99188460851979));
	rij_mins.insert(std::pair<std::string,double>("SER:SER",6.14628809264394));
	rij_mins.insert(std::pair<std::string,double>("SER:THR",6.40884644256299));
	rij_mins.insert(std::pair<std::string,double>("SER:TRP",8.12720237122348));
	rij_mins.insert(std::pair<std::string,double>("SER:TYR",8.13555345176175));
	rij_mins.insert(std::pair<std::string,double>("SER:VAL",6.60224140690369));
	rij_mins.insert(std::pair<std::string,double>("THR:THR",6.71212309439505));
	rij_mins.insert(std::pair<std::string,double>("THR:TRP",8.42519788658605));
	rij_mins.insert(std::pair<std::string,double>("THR:TYR",8.47061972948137));
	rij_mins.insert(std::pair<std::string,double>("THR:VAL",6.92600723741991));
	rij_mins.insert(std::pair<std::string,double>("TRP:TRP",9.96176991823581));
	rij_mins.insert(std::pair<std::string,double>("TRP:TYR",9.90096939046951));
	rij_mins.insert(std::pair<std::string,double>("TRP:VAL",8.80473773581833));
	rij_mins.insert(std::pair<std::string,double>("TYR:TYR",9.36930023419582));
	rij_mins.insert(std::pair<std::string,double>("TYR:VAL",8.40160233696928));
	rij_mins.insert(std::pair<std::string,double>("VAL:VAL",7.10756468404602));
}

double get_rij_min(const std::string &resnameA, const std::string &resnameB, double cutoff, std::map<std::string,double> &rij_mins){
	std::string key;
	key = resnameA+":"+resnameB;
	if (rij_mins.find (key) == rij_mins.end()){
		key = resnameB+":"+resnameA;
		if (rij_mins.find (key) == rij_mins.end()){
			return (cutoff);
		} 
		else {
			return (rij_mins.at(key));
		}
	} 
	else 
	{        
		return (rij_mins.at(key));
	}
}

int main (int argc, char **argv){
  int i;
  unsigned int f;
  std::stringstream resid;
  std::vector<std::string> pdbs;
  std::string currArg;
  std::string nucleus;
  std::string resname;
  std::string atomname;
  std::string frij_mins;
  std::string fcsfile;
  
  std::vector<std::string> trajs;
  std::vector<std::string> selections;
  std::vector<Molecule*> molecules;
  std::map<std::string,double> cs_data;
  std::map<std::string,std::string> nearset_neighbors;
  std::map<std::string,double> rij_mins;
  int start;
  int stop=std::numeric_limits<int>::max();
  int skip;
  bool startFlag=false;
  unsigned int itrj;
  std::ifstream trjin;
  Trajectory *ftrjin;
  unsigned int nframe;
  unsigned int process;
  double dist;
  double buffer;
  double cutoff;
  double rij_min;
  Molecule *Amol, *Bmol, *Cmol, *cmol;
  Atom *atmB, *atmC;  

  start=0;
  fcsfile="";
  frij_mins="";
  skip=0;
  nframe=0;
  process=0;
  f=0;    
  Amol=NULL;
  Bmol=NULL;
  Cmol=NULL;
  atmB=NULL;
  atmC=NULL;
  buffer=1.0;
  cutoff=10.0;
  rij_min=cutoff;
  pdbs.clear();
  selections.clear();
  molecules.clear();

  for (i=1; i<argc; i++){
    currArg=argv[i];
    if (currArg.compare("-h") == 0 || currArg.compare("-help") == 0)
    {
      usage();
    }
    else if (currArg.compare("-sel") == 0 )
    {
      currArg=argv[++i];    
      selections.push_back(currArg);
    }      
    else if (currArg.compare("-cutoff") == 0)
    {
      currArg=argv[++i];
      std::stringstream(currArg) >> cutoff;
    }
    else if (currArg.compare("-buffer") == 0)
    {
      currArg=argv[++i];
      std::stringstream(currArg) >> buffer;
    }
    else if (currArg.compare("-rij_mins") == 0 )
    {
      currArg=argv[++i];    
      frij_mins=currArg;     
    }          
    else if (currArg.compare("-csfile") == 0 )
    {
      currArg=argv[++i];    
      fcsfile=currArg;     
    }          
    else if (currArg.compare("-trj") == 0 || currArg.compare("-traj") == 0)
    {
      currArg=argv[++i];
      trajs.push_back(currArg);
    }
    else if (currArg.compare("-skip") == 0)
    {
      currArg=argv[++i];
      std::stringstream(currArg) >> skip;
    }
    else if (currArg.compare("-start") == 0)
    {
      currArg=argv[++i];
      std::stringstream(currArg) >> start;
      start--;
      startFlag=true;
    }
    else if (currArg.compare("-stop") == 0)
    {
      currArg=argv[++i];
      std::stringstream(currArg) >> stop;
    }    
    else if (currArg.compare(0,1,"-") == 0)
    {
      std::cerr << "Warning: Skipping unknown option \"" << currArg << "\"" << std::endl;
    }
    else{
      pdbs.push_back(currArg);
    }
  }
  if (pdbs.size() == 0)
  {
    std::cerr << std::endl << "Error: Please provide an input file" << std::endl << std::endl;
    usage();
  }
  if (pdbs.size() == 0)
  {
    std::cerr << std::endl << "Error: Please provide an input file" << std::endl << std::endl;
    usage();
  }
 
  if (selections.size() < 2 )
  {
    std::cerr << std::endl << " Error: supply at least two selection via -sel" << std::endl << std::endl;
    usage();
  }

  // loading rij_mins
  rij_mins.clear();
  if(frij_mins.length()==0){
  	load_rij_mins(rij_mins);
	}
	else
	{
		load_rij_mins_file(frij_mins,rij_mins);
	}

  // loading chemical shift data
  cs_data.clear();
  if(fcsfile.length()>0){
  	load_chem_shift_file(fcsfile,cs_data);
	}
	
	// generate nearest neighbors map
	cmol=Molecule::readPDB(pdbs.at(0));
	cmol->selAll();	
	nearset_neighbors.clear();
	setup_nearset_neighbors(cmol,nearset_neighbors);
  return(0);
  if (trajs.size() > 0)
  {
    if (pdbs.size() > 1)
    {
      std::cerr << std::endl  << "Warning: Only the first PDB structure is used for trajectory analysis" << std::endl << std::endl;
    }
    
    /* Trajectory analysis */
    cmol=Molecule::readPDB(pdbs.at(0));
    cmol->selAll();
        
    /* Process trajectories */
    for (itrj=0; itrj< trajs.size(); itrj++)
    {
      trjin.open(trajs.at(itrj).c_str(), std::ios::binary);
      if (trjin.is_open())
      {
        ftrjin=new Trajectory;
        ftrjin->setMolecule(cmol);
        if (ftrjin->findFormat(trjin) == true)
        {
          ftrjin->readHeader(trjin);
          if (skip > 0 && startFlag == false)
          {
            start=skip;
          }
          /* Loop through desired frames */
          for (i=start; i< ftrjin->getNFrame() && i< stop; i=i+1+skip)
          {
            if( ftrjin->readFrame(trjin, i) == false)
            {
              std::cerr << "Warning: EOF found before the next frame could be read" << std::endl;
              break;
            }
            nframe++;
            cmol->select(":.CA");
            Amol = cmol->copy(true);

						for (unsigned int k=0; k< selections.size(); k++)
						{
							Amol->select(selections.at(k));
							molecules.push_back(Amol->copy(true));
						}
                                  
            for (unsigned int j=0; j < molecules.size(); j++)
            {
							Bmol = molecules.at(j);
							for (unsigned int k=0; k < molecules.size(); k++)
							{								
								if(k>j)
								{
									Cmol = molecules.at(k);
									for (unsigned int l=0; l< Bmol->getAtmVecSize(); l++)
									{
										for (unsigned int m=0; m< Cmol->getAtmVecSize(); m++)
										{
											atmB = Bmol->getAtom(l);
											atmC = Cmol->getAtom(m);
											rij_min = buffer + get_rij_min(atmB->getResName(),atmC->getResName(),cutoff,rij_mins);
											dist = Analyze::distance(atmB->getCoor(), atmC->getCoor());					
											if (dist <= rij_min){ 
												 std::cout << nframe << " " << atmB->getResId() << " " << atmB->getResName() << " " << selections.at(j) << " " << atmC->getResId() << " "  << atmC->getResName() << " " << selections.at(k) << std::endl;
											}
										}
									}
								}
							}
            }
          }
        }
        else
        {
          std::cerr << "Warning: Skipping unknown trajectory format \"";
          std::cerr << trajs.at(itrj) << "\"" << std::endl;
        }
        if (ftrjin != NULL){
          delete ftrjin;
        }
      }
      trjin.close();
    }
  }
  else 
  { 
    for (f=0; f< pdbs.size(); f++)
    {  
      cmol = Molecule::readPDB(pdbs.at(f));
			cmol->select(":.CA");
			Amol = cmol->copy(true);

			for (unsigned int k=0; k< selections.size(); k++)
			{
				Amol->select(selections.at(k));
				molecules.push_back(Amol->copy(true));
				
			}
													
			for (unsigned int j=0; j < molecules.size(); j++)
			{
				Bmol = molecules.at(j);
				for (unsigned int k=0; k < molecules.size(); k++)
				{								
					if(k>j)
					{
						Cmol = molecules.at(k);
						for (unsigned int l=0; l< Bmol->getAtmVecSize(); l++)
						{
							for (unsigned int m=0; m< Cmol->getAtmVecSize(); m++)
							{
								atmB = Bmol->getAtom(l);
								atmC = Cmol->getAtom(m);
								rij_min = buffer + get_rij_min(atmB->getResName(),atmC->getResName(),cutoff,rij_mins);
								dist = Analyze::distance(atmB->getCoor(), atmC->getCoor());					
								if (dist <= rij_min){ 
									 std::cout << f+1 << " " << atmB->getResId() << " " << atmB->getResName() << " " << selections.at(j) << " " << atmC->getResId() << " "  << atmC->getResName() << " " << selections.at(k) << std::endl;
								}
							}
						}
					}
				}
			}
      delete cmol;
      delete Amol;
      molecules.clear();
    }
  }  
  return 0;
}
