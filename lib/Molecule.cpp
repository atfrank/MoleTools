//Sean M. Law

#include "Molecule.hpp"

int Molecule::readPDB (string *ifile, int model){

  ifstream pdbFile;
  string line;
  int currModel;
  double x,y,z;
  char lastChain;

  if (*ifile == "-"){

  }
  else{
    pdbFile.open((*ifile).c_str());
    if (pdbFile.is_open()){
      lastChain='0';
      while (pdbFile.good()){
        getline(pdbFile,line);
        //Check the model
        if (line.size() > 6 && line.compare(0,6,"MODEL ")==0){
          stringstream(line.substr(10,4)) >> currModel;
          if (model==0){
            model=1; //Use first model if undefined
          }
        }
        if (model && currModel != model){
          continue;
        }
        //Process atom entry
        if (line.size() >= 54 && (line.compare(0,4,"ATOM")==0 || line.compare(0,6,"HETATM")==0)){
          //substr: first character is denoted by a value of 0 (not 1)
          if ((line.substr(21,1))[0] != lastChain){
            lastChain=line.substr(21,1)[0];
            cout << "CHAIN " << lastChain << endl;
          }
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
          if (line.size() >= 60){
            stringstream(line.substr(54,6)) >> occu;
          }
          if (line.size() >= 66){
            stringstream(line.substr(60,6)) >> bfac;
          }
          if (line.size() >= 76){
            segid=line.substr(72,4);
          }
        */
          cout << x << " " << y << " " << z << endl;
        }
      }
      pdbFile.close();
    }
  }
  return 0;
}
