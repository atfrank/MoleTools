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
	int record;
	char header[4];

  trjin.seekg(0, std::ios::beg);
 
  trjin.read(reinterpret_cast<char*>(&record), sizeof(int));
	trjin.read(reinterpret_cast<char*>(&header), sizeof(char)*4);

	std::string hdr(header,4);

	if (record == 84 && hdr.compare("CORD") == 0){
		format="CHARMM";
    swab=false;
	}

  trjin.seekg(0, std::ios::beg);

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


void Trajectory::readHeader(std::ifstream &trjin){
	binbuf *buffer;
	int length;
	unsigned int i;

	buffer=NULL;

	trjin.seekg(0, std::ios::beg);

	if (format.compare("CHARMM") == 0){
		buffer=readFortran(trjin, buffer, length);

		//HDR
		for (i=0; i< 4; i++){
			hdr[i]=buffer[0].c[i];
		}

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

void Trajectory::readFrame(std::ifstream &trjin){
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
      for (i=0; i< natom; i++){
        mol->getAtom(i)->setCoor(Vector(xbuffer[i], ybuffer[i], zbuffer[i]));
      }
    }
    else{
		  for (i=0; i< natom; i++){
				x.push_back(xbuffer[i]);
			  y.push_back(ybuffer[i]);
			  z.push_back(zbuffer[i]);
		
			  std::cout << std::fixed;
			  std::cout << "coor" << std::setw(7) << i+1;
			  std::cout << std::setw(14) << xbuffer[i] << " ";
			  std::cout << std::setw(14) << ybuffer[i] << " ";
			  std::cout << std::setw(14) << zbuffer[i] << std::endl;
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
