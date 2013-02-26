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

Trajectory::binbuf* Trajectory::readFortran(std::ifstream &trjin){
	int recStart;
  int recEnd;
	binbuf *buffer;

	trjin.read(reinterpret_cast<char*>(&recStart), sizeof(int));
	buffer = new binbuf [recStart];
//	trjin.read((char*)(buffer), recStart);
	trjin.read(reinterpret_cast<char*>(buffer), recStart);
	trjin.read(reinterpret_cast<char*>(&recEnd), sizeof(int));

	return buffer;
}

void Trajectory::readHeader(std::ifstream &trjin){
	binbuf *bufout;
	unsigned int i;

	trjin.seekg(0, std::ios::beg);

	if (format.compare("CHARMM") == 0){
		bufout=readFortran(trjin);
		for (i=0; i< 4; i++){
			hdr[i]=bufout[0].c[i];
		}
		std::cerr << hdr << std::endl;
	}
	else if (format.compare("AMBER") == 0){
		
	}
	else{
		//Do nothing
	}
}

