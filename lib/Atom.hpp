//Sean M. Law

#include "Vector.hpp"

#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <cctype> //isdigit
#include <cstdlib>
using namespace std;

class Molecule {
  private:
    string recname; //Record name: "ATOM  ", "HETATM"
    int  atmnum; //Atom serial number
    string atmname; //Atom name
    char alt; //Alternate location indicator
    string resname; //Residue name
    char chainid; //Chain identifier
    int  resid; //Residue sequence number
    char icode; //Code for insertion of residues
    Vector coor; //X, Y, Z Coordinates
    double occu; //Occupancy
    double bfac; //B-factor or Temperature factor
    string segid; //Segment identifier
    int sel; //Selection flag

  public:
    Molecule(); //Constructor
    Molecule(int atmnum, string atmname, string resname, int resnum, Vector vec, string seg=0); //Overload constructor

    int readPDB (string *ifile);
};
