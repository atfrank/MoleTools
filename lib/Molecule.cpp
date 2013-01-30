//Sean M. Law

#include "Molecule.hpp"

int Atom::readPDB (string *ifile){

  ifstream pdbFile;
  string line;
  double x,y,z;

  if (*ifile == "-"){

  }
  else{
    pdbFile.open((*ifile).c_str());
    if (pdbFile.is_open()){
      while (pdbFile.good()){
        getline(pdbFile,line);
        if (line.length() > 50 && line.compare(1,4,"ATOM") && line.compare(1,6,"HETATM")){
/*
          stringstream(line.substr(7,5)) >> atmnum;
          atmname=line.substr(13,4);
          alt=(line.substr(17,1))[0];
          resname=line.substr(18,3);
          chainid=(line.substr(22,1))[0];
          stringstream(line.substr(23,4)) >> resid;
          icode=(line.substr(27,1))[0];
          stringstream(line.substr(31,8)) >> x;
          stringstream(line.substr(39,8)) >> y;
          stringstream(line.substr(47,8)) >> z;
          coor=Vector(x,y,z);
          stringstream(line.substr(55,6)) >> occu;
          stringstream(line.substr(61,6)) >> bfac;
          segid=line.substr(73,4);
*/
        }
      }
      pdbFile.close();
    }
  }
  return 0;
}
