//Sean M. Law

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
    char alt; //Alternate location indicator
    std::string resname; //Residue name
    char chainid; //Chain identifier
    int  resid; //Residue sequence number
    char icode; //Code for insertion of residues
    Vector coor; //X, Y, Z Coordinates
    double occu; //Occupancy
    double bfac; //B-factor or Temperature factor
    std::string segid; //Segment identifier
    int sel; //Selection flag

  public:
    Atom(); //Constructor
    Atom(int atmnum, std::string atmname, std::string resname, int resnum, Vector vec, std::string seg=0); //Overload constructor

    //Get atom info
    std::string& getRecName();
    int& getAtmNum();
    std::string& getAtmName();
    char& getAlt();
    std::string& getResName();
    char& getChainId();
    int& getResId();
    char& getICode();
    Vector& getCoor ();
    double& getX();
    double& getY();
    double& getZ();
    double& getOccu();
    double& getBFac();
    std::string& getSegId();

    //Set atom info
    void setRecName(const std::string& recnamein);
    void setAtmNum(const int& atmnumin);
    void setAtmName(const std::string& atmnamein);
    void setAlt(const char& altin);
    void setResName(const std::string& resnamein);
    void setChainId(const char& chainidin);
    void setResId(const int& residin);
    void setICode(const char& icodein);
    void setCoor(const Vector& coorin);
    void setOccu(const double& occuin);
    void setBFac(const double& bfacin);
    void setSegId(const std::string& segidin);
};
