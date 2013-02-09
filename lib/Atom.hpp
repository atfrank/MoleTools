//Sean M. Law
#ifndef ATOM_H
#define ATOM_H

#include "Vector.hpp"

#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <cctype> //isdigit
#include <cstdlib>

class Atom {
  private:
    std::string recname; //Record name: "ATOM  ", "HETATM"
    int  atmnum; //Atom serial number
    std::string atmname; //Atom name
    std::string alt; //Alternate location indicator
    std::string resname; //Residue name
    std::string chainid; //Chain identifier, modified if realid is duplicated
    std::string realid; //Store original chainid
    int  resid; //Residue sequence number
    std::string icode; //Code for insertion of residues
    Vector coor; //X, Y, Z Coordinates
    double occu; //Occupancy
    double bfac; //B-factor or Temperature factor
    std::string segid; //Segment identifier
    bool sel; //Selection flag

  public:
    Atom(); //Constructor
    Atom(int atmnum, std::string atmname, std::string resname, int resnum, Vector vec, std::string seg=0); //Overload constructor

    void reset();
    void clone(Atom* atmin);

    //Get atom info
    std::string& getRecName();
    int& getAtmNum();
    std::string& getAtmName();
    std::string& getAlt();
    std::string& getResName();
    std::string& getChainId();
    std::string& getRealId();
    int& getResId();
    std::string& getICode();
    Vector& getCoor ();
    double& getX();
    double& getY();
    double& getZ();
    double& getOccu();
    double& getBFac();
    std::string& getSegId();
    bool& getSel();

    //Set atom info
    void setRecName(const std::string& recnamein);
    void setAtmNum(const int& atmnumin);
    void setAtmName(const std::string& atmnamein);
    void setAtmName(); //Clear
    void setAlt(const std::string& altin);
    void setAlt(); //Clear
    void setResName(const std::string& resnamein);
    void setResName(); //Clear
    void setChainId(const std::string& chainidin);
    void setRealId(const std::string& realidin);
    void setChainId(); //Clear
    void setRealId(); //Clear
    void setResId(const int& residin);
    void setICode(const std::string& icodein);
    void setICode(); //Clear
    void setCoor(const Vector& coorin);
    void setOccu(const double& occuin);
    void setBFac(const double& bfacin);
    void setSegId(const std::string& segidin);
    void setSegId(); //Clear
    void setSel(const bool selin);
};

#endif
