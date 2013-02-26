//Sean M. Law

#include "Trajectory.hpp"

//Deal with Endianess

bool Trajectory::findFormat(std::ifstream &trjin){
	int record;
	char header[4];

  trjin.seekg(0, std::ios::beg);
 
  trjin.read(reinterpret_cast<char*>(&record), sizeof(int));
	trjin.read(reinterpret_cast<char*>(&header), sizeof(char)*4);

	std::string hdr(header);

	if (record == 84 && hdr.compare("CORD") == 0){
		format="CHARMM";
	}

  trjin.seekg(0, std::ios::beg);

	if (format.length() > 0){
		return true;
	}
	else{
		return false;
	}
}

Trajectory::binbuf* Trajectory::readFortran(std::ifstream &trjin, int &length){
	int recStart;
  int recEnd;
	binbuf *buffer;

	//Read Fortran record lengths and buffer
	trjin.read(reinterpret_cast<char*>(&recStart), sizeof(int));
	buffer = new binbuf [recStart];
	trjin.read(reinterpret_cast<char*>(buffer), recStart);
	trjin.read(reinterpret_cast<char*>(&recEnd), sizeof(int));

	
	//Check Fortran record length mismatch
	if (recStart == recEnd){
		length=recStart;
	}
	else{
		std::cerr << "Error: Fortran record length marker mismatch" << std::endl;
		std::exit(0);
	}

	return buffer;
}

void Trajectory::readHeader(std::ifstream &trjin){
	binbuf *buffer;
	int length;
	unsigned int i;

	buffer=NULL;

	trjin.seekg(0, std::ios::beg);

	if (format.compare("CHARMM") == 0){
		buffer=readFortran(trjin, length);
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
		else{
			fixed=false;
		}
		tstep=AKMATPS*buffer[10].f;
		
		if (buffer[11].i > 0){
			crystal=true;
		}
		else{
			crystal=false;
		}	
		
		//first=tstart/delta; //Divided by zero!!
		deltat=delta*tstep;

		//Title
		buffer=readFortran(trjin, length);

		//NATOM
		buffer=readFortran(trjin, length);
		natom=buffer[0].i;

    //FIXED
    if (fixed == true){
      buffer=readFortran(trjin, length);
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
	binbuf *buffer;
	binbuf *xbuffer;
	binbuf *ybuffer;
	binbuf *zbuffer;
  int length;
	int i;

	buffer=NULL;
	xbuffer=NULL;
	ybuffer=NULL;
	zbuffer=NULL;

	if (crystal == true){
		buffer=readFortran(trjin, length);
		if (buffer != NULL){
			delete buffer;
		}
	}

	//Coordinates
	xbuffer=readFortran(trjin, length);
	ybuffer=readFortran(trjin, length);
	zbuffer=readFortran(trjin, length);

	x.clear();
	y.clear();
	z.clear();

	if (fixed == true){
		std::cerr << "Warning: Fixed atoms has yet to be implemented"<< std::endl;
	}
	else{
		for (i=0; i< natom; i++){
			x.push_back(xbuffer[i].f);
			y.push_back(ybuffer[i].f);
			z.push_back(zbuffer[i].f);
			//std::cout << xbuffer[i].f << "  " << ybuffer[i].f << "  " << zbuffer[i].f << std::endl;
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
