//Sean M. Law

#include <iomanip>

#include "Atom.hpp"

Atom::Atom(){
  recname="ATOM";
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

Atom::Atom(int atmnumin, std::string atmnamein, std::string resnamein, int residin, Vector coorin, std::string segidin){

  recname="ATOM";
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

//Get atom info
std::string& Atom::getRecName(){
  return recname;
}

int& Atom::getAtmNum(){
  return atmnum;
}

std::string& Atom::getAtmName(){
  return atmname;
}

char& Atom::getAlt(){
  return alt;
}

std::string& Atom::getResName(){
  return resname;
}

char& Atom::getChainId(){
  return chainid;
}

int& Atom::getResId(){
  return resid;
}

char& Atom::getICode(){
  return icode;
}

Vector& Atom::getCoor () {
  return coor;
}

double& Atom::getX () {
  return coor.x();
}

double& Atom::getY () {
  return coor.y();
}

double& Atom::getZ () {
  return coor.z();
}

double& Atom::getOccu(){
  return occu;
}

double& Atom::getBFac(){
  return bfac;
}

std::string& Atom::getSegId(){
  return segid;
}

//Set atom info
void Atom::setRecName(const std::string& recnamein){
  this->recname=recnamein;
}

void Atom::setAtmNum(const int& atmnumin){
  this->atmnum=atmnumin;
}

void Atom::setAtmName(const std::string& atmnamein){
  this->atmname=atmnamein;
}

void Atom::setAlt(const char& altin){
  this->alt=altin;
}

void Atom::setResName(const std::string& resnamein){
  this->resname=resnamein;
}

void Atom::setChainId(const char& chainidin){
  this->chainid=chainidin;
}

void Atom::setResId(const int& residin){
  this->resid=residin;
}

void Atom::setICode(const char& icodein){
  this->icode=icodein;
}

void Atom::setCoor (const Vector& coorin){
  this->coor=coorin;
}

void Atom::setOccu(const double& occuin){
  this->occu=occuin;
}

void Atom::setBFac(const double& bfacin){
  this->bfac=bfacin;
}

void Atom::setSegId(const std::string& segidin){
  this->segid=segidin;
}

