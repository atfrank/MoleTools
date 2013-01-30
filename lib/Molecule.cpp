//Sean M. Law

#include "Molecule.hpp"

int Molecule::readPDB (string *ifile){

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
        if (line.size() > 50 && (line.compare(0,4,"ATOM")==0 || line.compare(0,6,"HETATM")==0)){
          //substr: first character is denoted by a value of 0 (not 1)
        /*
          stringstream(line.substr(6,5)) >> atmnum;
          atmname=line.substr(12,4);
          alt=(line.substr(16,1))[0];
          resname=line.substr(17,3);
          chainid=(line.substr(21,1))[0];
          stringstream(line.substr(22,4)) >> resid;
          icode=(line.substr(26,1))[0];
        */
          stringstream(line.substr(30,8)) >> x;
          stringstream(line.substr(38,8)) >> y;
          stringstream(line.substr(46,8)) >> z;
        /*
          coor=Vector(x,y,z);
          stringstream(line.substr(54,6)) >> occu;
          stringstream(line.substr(60,6)) >> bfac;
          segid=line.substr(72,4);
        */
          cout << x << " " << y << " " << z << endl;
        }
      }
      pdbFile.close();
    }
  }
  return 0;
}
