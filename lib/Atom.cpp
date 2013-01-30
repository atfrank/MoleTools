#include "Atom.hpp"

Atom::Atom(){
  atmnum=0;
  atmname="0";
  alt=0;
  resname="0";
  chainid=0;
  resid=0;
  coor=Vector(0.0, 0.0, 0.0);
  occu=0.0;
  bfac=0.0;
  segid="0";
  sel=0;
}

Atom::Atom(int atmnumin, string atmnamein, string resnamein, int residin, Vector coorin, string segidin){

  atmnum=atmnumin;
  atmname=atmnamein;
  alt=0;
  resname=resnamein;
  chainid=0;
  resid=residin;
  coor=coorin;
  occu=0.0;
  bfac=0.0;
  segid=segidin;
  sel=0;
}

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
        }
      }
      pdbFile.close();
    }
  }
  return 0;
}
