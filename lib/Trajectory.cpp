//Sean M. Law

#include "Trajectory.hpp"

//Deal with Endianess

Trajectory::Trajectory (){
  swab=false;
  mol=NULL;
  hdr.clear();
  nframe=0;
  npriv=0;
  nsavc=0;
  nstep=0;
  qvelocity=false;
  dof=0;
  nfixed=0;
  tstep=0.0;
  qcrystal=false;
  q4d=false;
  qcharge=false;
  qcheck=false;
  version=0;
  title1.clear();
  title2.clear();
  natom=0;
  fixinx.clear();
}

bool Trajectory::findFormat(std::ifstream &trjin){
 	binbuf *buffer;
  int length;

  buffer=NULL;

  trjin.seekg(0, std::ios::beg);

  buffer=readFortran(trjin, buffer, length);

	hdr.assign(buffer[0].c,4);

	if (length == 84 && hdr.compare("CORD") == 0){
		format="CHARMM";
    swab=false;
	}

  trjin.seekg(0, std::ios::beg);

	clearHeader();

	if (format.length() > 0){
		return true;
	}
	else{
		return false;
	}
}

template <class BinBuf> 
BinBuf* Trajectory::readFortran(std::ifstream &trjin, BinBuf *buffer, int &length){
  int recStart;
  int recEnd;
	BinBuf *binOut;

  //Read Fortran record lengths and buffer
  trjin.read(reinterpret_cast<char*>(&recStart), sizeof(int));
  binOut = new BinBuf [recStart];
  trjin.read(reinterpret_cast<char*>(binOut), recStart);
  trjin.read(reinterpret_cast<char*>(&recEnd), sizeof(int));

  //Check Fortran record length mismatch
  if (recStart == recEnd){
    length=recStart;
  }
  else{
    std::cerr << "Error: Fortran record length marker mismatch" << std::endl;
    std::exit(0);
  }

	return binOut;
}

template <class BinBuf>
void Trajectory::writeFortran(std::ofstream &trjout, BinBuf *buffer, int &length){
  trjout.seekp(0, std::ios::beg);
  //Write Fortran record lengths and buffer
  trjout.write(reinterpret_cast<char*>(&length), sizeof(int));
  //binOut = new BinBuf [recStart];
  trjout.write(reinterpret_cast<char*>(buffer), length);
  trjout.write(reinterpret_cast<char*>(&length), sizeof(int));
}

void Trajectory::clearHeader(){
	hdr.clear();
	nframe=0;
	npriv=0;
	nsavc=0;
	nstep=0;
	qvelocity=0;
	dof=0;
	nfixed=0;
	tstep=0.0;
	qcrystal=0;
	q4d=0;
	qcharge=0;
	qcheck=0;
	version=0;
	title1.clear();
  title2.clear();
	natom=0;
	fixinx.clear();
}

void Trajectory::readHeader(std::ifstream &trjin){
	binbuf *buffer;
	char *cbuffer;
	int length;
	//int i;

	buffer=NULL;
	cbuffer=NULL;

	trjin.seekg(0, std::ios::beg);

	if (format.compare("CHARMM") == 0){
		buffer=readFortran(trjin, buffer, length);

		//HDR
		hdr.assign(buffer[0].c,4);

		//ICNTRL 1-20
		nframe=buffer[1].i;
		npriv=buffer[2].i;
		nsavc=buffer[3].i;
		nstep=buffer[4].i;
		qvelocity=buffer[5].i;
    if (qvelocity > 0){
			std::cerr << "Error: Velocity reading has yet to be implemented" << std::endl;
		}
	
		dof=buffer[8].i;
		nfixed=buffer[9].i;
	
		tstep=buffer[10].f;
		
		qcrystal=buffer[11].i;
		q4d=buffer[12].i;
		qcharge=buffer[13].i;
		qcheck=buffer[14].i;
		
		version=buffer[20].i;
		
		//Title
		cbuffer=readFortran(trjin, cbuffer, length);
		title1.assign(cbuffer,80);
    //81-84 are blank spaces, NOT null characters (\0)
		title2.assign(cbuffer+84,80);
		

		//NATOM
		buffer=readFortran(trjin, buffer, length);
		natom=buffer[0].i;

    //FIXED
    if (nfixed > 0){
      buffer=readFortran(trjin, buffer, length);
      std::cerr << "Warning: Fixed atoms has yet to be implemented" << std::endl;
      //for (i=0; i< length; i++){
			//	std::cerr << buffer[i].i << ":";
      //  fixinx.push_back(buffer[i].i);
      //}
    }
	
	}
	else if (format.compare("AMBER") == 0){
		
	}
	else{
		//Do nothing
	}

	if (buffer != NULL){
		delete buffer;
	}
	if (cbuffer != NULL){
		delete cbuffer;
	}
}

void Trajectory::writeHeader(std::ofstream &trjout){
  if (format == "CHARMM"){
    unsigned int icntrl[21];
    int length;
    binbuf *buffer=reinterpret_cast<binbuf *>(&icntrl[0]);
    for (unsigned int i=0; i< 21; i++){
      icntrl[i]=0;
    }

    buffer[0].c[0]='C';
    buffer[0].c[1]='O';
    buffer[0].c[2]='R';
    buffer[0].c[3]='D';

    icntrl[1]=getNFrame();
    icntrl[2]=getNPriv();
    icntrl[3]=getNSavc();
    icntrl[4]=getNStep();
    icntrl[5]=getQVelocity();

    icntrl[8]=getDOF();
    icntrl[9]=getNFixed();
    buffer[10].f=getTStepAKMA();
    icntrl[11]=getQCrystal();
    icntrl[12]=getQ4D();
    icntrl[13]=getQCharge();
    icntrl[14]=getQCheck();

    icntrl[20]=getVersion();

    length=sizeof(int)*21;

    writeFortran(trjout, buffer, length);
  }
}

void Trajectory::showHeader(){
	std::cerr << title1 << std::endl;
  std::cerr << title2 << std::endl;
	std::cerr << std::fixed;
  std::cerr << std::setw(25) << std::left << "Atoms" << ": " << natom << std::endl;
	std::cerr << std::setw(25) << std::left << "Frames" << ": " << nframe << std::endl;
  std::cerr << std::setw(25) << std::left << "Start Frame" << ": " << npriv << std::endl;
  std::cerr << std::setw(25) << std::left << "Save Frequency" << ": " << nsavc << std::endl;
  std::cerr << std::setw(25) << std::left << "Dynamics Steps" << ": " << nstep << std::endl;
  std::cerr << std::setw(25) << std::left << "Degrees of Freedom" << ": " << dof << std::endl;
	std::cerr << std::setw(25) << std::left << "Number of Fixed" << ": " << nfixed << std::endl;
  std::cerr << std::setw(25) << std::left << "Time Step (ps)" << ": " << this->getTStepPS() << std::endl;
	std::cerr << std::setw(25) << std::left << "Start Time (ps)" << ": " << (npriv/nsavc)*this->getTStepPS() <<  std::endl;
  std::cerr << std::setw(25) << std::left << "Periodic Boundaries" << ": " << qcrystal << std::endl;
  std::cerr << std::setw(25) << std::left << "4D Trajectory" << ": " << q4d << std::endl;
  std::cerr << std::setw(25) << std::left << "Fluctuating Charges" << ": " << qcharge << std::endl;
  std::cerr << std::setw(25) << std::left << "Consistency Check" << ": " << qcheck << std::endl;
  std:: cerr << std::setw(25) << std::left << "Version" << ": " << version << std::endl;
}

void Trajectory::cloneHeader(Trajectory *ftrjin){
  format=ftrjin->getFormat();
	hdr=ftrjin->getHdr();
  nframe=ftrjin->getNFrame();
  npriv=ftrjin->getNPriv();
  nsavc=ftrjin->getNSavc();
  nstep=ftrjin->getNStep();
  qvelocity=ftrjin->getQVelocity();

  dof=ftrjin->getDOF();
  nfixed=ftrjin->getNFixed();
  tstep=ftrjin->getTStepAKMA();
  qcrystal=ftrjin->getQCrystal();
  q4d=ftrjin->getQ4D();
  qcharge=ftrjin->getQCharge();
  qcheck=ftrjin->getQCheck();

  version=ftrjin->getVersion();

  title1=ftrjin->getTitle1();
  title2=ftrjin->getTitle2();
  natom=ftrjin->getNAtom();
	for (unsigned int i=0; i< ftrjin->getFixInxVecSize(); i++){
   	fixinx.push_back(ftrjin->getFixInx(i));
	}
}

void Trajectory::readFrame(std::ifstream &trjin, unsigned int frame){
	double *dbuffer; //Needed for crystal!
	float *xbuffer;
	float *ybuffer;
	float *zbuffer;
  int length;
	int i;

	dbuffer=NULL;
	xbuffer=NULL;
	ybuffer=NULL;
	zbuffer=NULL;

	if (qcrystal == true){
		dbuffer=readFortran(trjin, dbuffer, length);

		//Box dimensions
		pbx=sqrt(dbuffer[0]*dbuffer[0]+dbuffer[1]*dbuffer[1]+dbuffer[3]*dbuffer[3]);
    pby=sqrt(dbuffer[1]*dbuffer[1]+dbuffer[2]*dbuffer[2]+dbuffer[4]*dbuffer[4]);
    pbz=sqrt(dbuffer[3]*dbuffer[3]+dbuffer[4]*dbuffer[4]+dbuffer[5]*dbuffer[5]);

		//Box angles
    pbalpha=acos(dbuffer[4]*(dbuffer[2]+dbuffer[5])+dbuffer[1]*dbuffer[3]/(pby*pbz))*180.0/PI;
    pbbeta=acos(dbuffer[3]*(dbuffer[0]+dbuffer[5])+dbuffer[1]*dbuffer[4]/(pbz*pbx))*180.0/PI;
    pbgamma=acos(dbuffer[1]*(dbuffer[0]+dbuffer[2])+dbuffer[3]*dbuffer[4]/(pbx*pby))*180.0/PI;
		
		/*
		std::cerr << "Periodic Box :";
		std::cerr << pbx << " x " << pby << " x " << pbz << " ";
		std::cerr << "( " << pbalpha << " x " << pbbeta << " x " << pbgamma << " )" << std::endl;
		*/
		
		delete dbuffer;
	}

	//Coordinates
	xbuffer=readFortran(trjin, xbuffer, length);
	ybuffer=readFortran(trjin, ybuffer, length);
	zbuffer=readFortran(trjin, zbuffer, length);

	if (nfixed > 0){
		std::cerr << "Warning: Fixed atoms has yet to be implemented"<< std::endl;
	}
	else{
    if (mol != NULL){
			if (length/sizeof(float) != mol->getNAtom()){
				std::cerr << "Error: The natom mismatch!" << std::endl;
			}
			else{
      	for (i=0; i< natom; i++){
        	mol->getAtom(i)->setCoor(Vector(static_cast<double>(xbuffer[i]), static_cast<double>(ybuffer[i]), static_cast<double>(zbuffer[i])));
      	}
			}
    }
    else{
		  for (i=0; i< natom; i++){
				x.push_back(xbuffer[i]);
			  y.push_back(ybuffer[i]);
			  z.push_back(zbuffer[i]);
			
			/*	
			  std::cout << std::fixed;
			  std::cout << "coor" << std::setw(7) << i+1;
			  std::cout << std::setw(14) << xbuffer[i] << " ";
			  std::cout << std::setw(14) << ybuffer[i] << " ";
			  std::cout << std::setw(14) << zbuffer[i] << std::endl;
			*/
      }
		}
	}

	if (xbuffer != NULL){
		delete xbuffer;
	}
	if (ybuffer != NULL){
		delete ybuffer;
	}
	if (zbuffer != NULL){
		delete zbuffer;
	}
}

void Trajectory::setMolecule(Molecule *molin){
  mol=molin;
}

std::string Trajectory::getFormat(){
  return format;
}

std::string Trajectory::getHdr(){
	return hdr;
}

int Trajectory::getNFrame(){
	return nframe;
}

int Trajectory::getNPriv(){
	return npriv;
}

int Trajectory::getNSavc(){
	return nsavc;
}

int Trajectory::getNStep(){
	return nstep;
}

int Trajectory::getQVelocity(){
	return qvelocity;
}

int Trajectory::getDOF(){
	return dof;
}
int Trajectory::getNFixed(){
	return nfixed;
}

float Trajectory::getTStepAKMA(){
	return tstep;
}

double Trajectory::getTStepPS(){
  return  static_cast<double>(static_cast<int>(static_cast<double>(tstep)*AKMATPS*nsavc*100000.0+0.5))/100000.0;
}

int Trajectory::getQCrystal(){
	return qcrystal;
}

int Trajectory::getQ4D(){
	return q4d;
}

int Trajectory::getQCharge(){
	return qcharge;
}

int Trajectory::getQCheck(){
	return qcheck;
}

int Trajectory::getVersion(){
	return version;
}

std::string Trajectory::getTitle1(){
  if (title1.find_first_of("*") != 0){
    title1.insert(0,"*");
  }
  title1.resize(80, ' ');
	return title1;
}

std::string Trajectory::getTitle2(){
  if (title2.find_first_of("*") != 0){
    title2.insert(0,"*");
  }
  title2.resize(80, ' ');
  return title2;
}

int Trajectory::getNAtom(){
	return natom;
}

unsigned int Trajectory::getFixInxVecSize(){
	return fixinx.size();
}

int Trajectory::getFixInx(int element){
	return fixinx.at(element);
}

void Trajectory::setHdr(const std::string &hdrin){
	hdr=hdrin;
}
void Trajectory::setNFrame(const int &nframein){
	nframe=nframein;
}

void Trajectory::setNPriv(const int &nprivin){
	npriv=nprivin;
}

void Trajectory::setNSavc(const int &nsavcin){
	nsavc=nsavcin;
}

void Trajectory::setNStep(const int &nstepin){
	nstep=nstepin;
}

void Trajectory::setQVelocity(const int &qvelocityin){
	qvelocity=qvelocityin;
}

void Trajectory::setDOF(const int &dofin){
	dof=dofin;
}

void Trajectory::setNFixed(const int &nfixedin){
	nfixed=nfixedin;
}

void Trajectory::setTStep(const double &tstepin){
  //Convert from picoseconds double to AKMA float
  tstep=static_cast<float>((static_cast<double>(static_cast<int>(tstepin*100000.0))-0.5)/100000.0/nsavc/AKMATPS);
}

void Trajectory::setTStep(const float &tstepin){
  //Float input is already in AKMA units
  tstep=tstepin;
}

void Trajectory::setQCrystal(const int &qcrystalin){
	qcrystal=qcrystalin;
}

void Trajectory::setQ4D(const int &q4din){
	q4d=q4din;
}

void Trajectory::setQCharge(const int &qchargein){
	qcharge=qchargein;
}

void Trajectory::setQCheck(const int &qcheckin){
	qcheck=qcheckin;
}

void Trajectory::setVersion(const int &versionin){
	version=versionin;
}

void Trajectory::setTitle1(const std::string &title1in){
	title1=title1in;
  if (title1.find_first_of("*") != 0){
    title1.insert(0,"*");
  }
  title1.resize(80, ' ');
}

void Trajectory::setTitle2(const std::string &title2in){
  title2=title2in;
  if(title2.find_first_of("*") != 0){
    title2.insert(0,"*");
  }
  title2.resize(80, ' ');
}

void Trajectory::setNAtom(const int &natomin){
	natom=natomin;
}
void Trajectory::addFixInx(const int &elementin){
	fixinx.push_back(elementin);
}

void Trajectory::clearFixInx(){
	fixinx.clear();
}
