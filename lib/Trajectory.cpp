//Sean M. Law

#include "Trajectory.hpp"

//Deal with Endianess

Trajectory::Trajectory (){
  swab=false;
  mol=NULL;
  crystal=false;
  fixed=false;
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
	tstart=0;
	first=0;
	delta=0;
	deltat=0.0;
	dof=0;
	tstep=0.0;
	crystal=false;
	pbx=0.0;
	pby=0.0;
	pbz=0.0;
	fixed=false;
	version=0;
	natom=0;
	endian.clear();
	title.clear();
	fixinx.clear();
}

void Trajectory::readHeader(std::ifstream &trjin){
	binbuf *buffer;
	int length;
	unsigned int i;

	buffer=NULL;

	trjin.seekg(0, std::ios::beg);

	if (format.compare("CHARMM") == 0){
		buffer=readFortran(trjin, buffer, length);

		//HDR
		hdr.assign(buffer[0].c,4);

		//ICNTRL 1-20
		nframe=buffer[1].i;
		tstart=buffer[2].i;
		delta=buffer[3].i;
		dof=buffer[8].i;
		if (buffer[9].i > 0){
			fixed=true;
		}
	
		tstep=AKMATPS*buffer[10].f;
		
		if (buffer[11].i > 0){
			crystal=true;
		}
		
		//first=tstart/delta; //Divided by zero!!
		deltat=delta*tstep;

		//Title
		buffer=readFortran(trjin, buffer, length);
		for (i=0; i< static_cast<unsigned int>(length); i++){
			
		}

		//NATOM
		buffer=readFortran(trjin, buffer, length);
		natom=buffer[0].i;

    //FIXED
    if (fixed == true){
      buffer=readFortran(trjin, buffer, length);
      std::cerr << "Warning: Fixed atoms has yet to be implemented"<< std::endl;
      //for (i=0; i< length; i++){
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
}

void Trajectory::readFrame(std::ifstream &trjin, unsigned int frame){
	binbuf *buffer; //Needed for crystal!
	float *xbuffer;
	float *ybuffer;
	float *zbuffer;
  int length;
	int i;

	buffer=NULL;
	xbuffer=NULL;
	ybuffer=NULL;
	zbuffer=NULL;

	if (crystal == true){
		buffer=readFortran(trjin, buffer, length);
		if (buffer != NULL){
			delete buffer;
		}
	}

	//Coordinates
	xbuffer=readFortran(trjin, xbuffer, length);
	ybuffer=readFortran(trjin, ybuffer, length);
	zbuffer=readFortran(trjin, zbuffer, length);

	if (fixed == true){
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
