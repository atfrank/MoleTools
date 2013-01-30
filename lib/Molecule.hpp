//Sean M. Law

#include "Vector.hpp"

#include <fstream>
#include <iostream>
#include <string>
#include <cstdlib>
using namespace std;

class Molecule {
  private:
    char recname[6]; //Record name: "ATOM  ", "HETATM"
    int  atmnum; //Atom serial number
    char atmname[6]; //Atom name
    char alt; //Alternate location indicator
    char resname[6]; //Residue name
    char chainid; //Chain identifier
    int  resid; //Residue sequence number
    char icode; //Code for insertion of residues
    Vector coor; //X, Y, Z Coordinates
    double occu; //Occupancy
    double bfac; //B-factor or Temperature factor
    char segid[6]; //Segment identifier
    int sel; //Selection flag

  public:
    Molecule(); //Constructor
    Molecule(int atmnum, char *atmname, char *resname, int resnum, Vector vec, char *seg=0); //Overload constructor

    void readPDB (string *ifile);
};
