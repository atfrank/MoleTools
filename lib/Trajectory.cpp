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
  title.clear();
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

void Trajectory::clearHeader(){
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
	title.clear();
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
		if (buffer[5].i > 0){
			qvelocity=true;
			std::cerr << "Error: Velocity reading has yet to be implemented" << std::endl;
		}
		else{
			qvelocity=false;
		}
	
		dof=buffer[8].i;
		nfixed=buffer[9].i;
	
		tstep=static_cast<double>(buffer[10].f);
		
		if (buffer[11].i > 0){
			qcrystal=true;
		}
		else{
			qcrystal=false;
		}
		if (buffer[12].i > 0){
			q4d=true;
		}
		else{
			q4d=false;
		}
		if (buffer[13].i > 0){
			qcharge=true;
		}
		else{
			qcharge=false;
		}
		if (buffer[14].i > 0){
			qcheck=true;
		}
		else{
			qcheck=false;
		}
		
		version=buffer[20].i;
		
		//Title
		cbuffer=readFortran(trjin, cbuffer, length);
		title.assign(cbuffer,80);
		title.append("\n");
		title.append(cbuffer+84,80);
		

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

void Trajectory::showHeader(){
	std::cerr << title << std::endl;
	std::cerr << std::fixed;
	std::cerr << std::setw(25) << std::left << "Atoms" << ": " << natom << std::endl;
	std::cerr << std::setw(25) << std::left << "Frames" << ": " << nframe << std::endl;
  std::cerr << std::setw(25) << std::left << "Degrees of Freedom" << ": " << dof << std::endl;
	std::cerr << std::setw(25) << std::left << "Fixed" << ": " << nfixed << std::endl;
  std::cerr << std::setw(25) << std::left << "Time Step" << ": " << static_cast<double>(static_cast<int>(tstep*AKMATPS*nsavc*100000.0+0.5))/100000.0 << std::endl;
	std::cerr << std::setw(25) << std::left << "Start Time" << ": " << (npriv/nsavc)*(static_cast<double>(static_cast<int>(tstep*AKMATPS*nsavc*100000.0+0.5))/100000.0) << std::endl;
  std::cerr << std::setw(25) << std::left << "Periodic Boundaries" << ": " << qcrystal << std::endl;
  std::cerr << std::setw(25) << std::left << "4D Trajectory" << ": " << q4d << std::endl;
  std::cerr << std::setw(25) << std::left << "Fluctuating Charges" << ": " << qcharge << std::endl;
  std::cerr << std::setw(25) << std::left << "Consistency Check" << ": " << qcheck << std::endl;
  std:: cerr << std::setw(25) << std::left << "Version" << ": " << version << std::endl;
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
			if (length/sizeof(float) != mol->getNatom()){
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

int Trajectory::getNFrame(){
	return nframe;
}

void Trajectory::setMolecule(Molecule *molin){
  mol=molin;
}
