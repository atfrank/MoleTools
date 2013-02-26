//Sean M. Law

#include "Trajectory.hpp"

//Deal with Endianess

void Trajectory::getFormat(std::ifstream &traj){
  char *buffer;
  int size;

  traj.seekg(0, std::ios::beg);
  
  size=1;
  buffer=new char [size];
  traj.read(buffer, size);



  traj.seekg(0, std::ios::beg);
}

std::string Trajectory::readFortran(){

}
