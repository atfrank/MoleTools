//Sean M. Law

#include <iomanip>

#include "Atom.hpp"

Atom::Atom(){
  recname="ATOM";
  atmnum=0;
  atmname.clear();
  alt.clear();
  resname.clear();
  chainid.clear();
  resid=0;
  coor=Vector(0.0, 0.0, 0.0);
  occu=0.0;
  bfac=0.0;
  segid.clear();
  sel=0;
}

Atom::Atom(int atmnumin, std::string atmnamein, std::string resnamein, int residin, Vector coorin, std::string segidin){

  recname="ATOM";
  atmnum=atmnumin;
  atmname=atmnamein;
  alt.clear();
  resname=resnamein;
  chainid.clear();
  resid=residin;
  coor=coorin;
  occu=0.0;
  bfac=0.0;
  segid=segidin;
  sel=0;
}

void Atom::reset(){
  this->setRecName("ATOM");
  this->setAtmNum(0);
  this->setAtmName();
  this->setAlt();
  this->setResName();
  this->setChainId();
  this->setResId(0);
  this->setCoor(Vector(0.0, 0.0, 0.0));
  this->setOccu(0.0);
  this->setBFac(0.0);
  this->setSegId();
  this->setSel(0);  
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

std::string& Atom::getAlt(){
  return alt;
}

std::string& Atom::getResName(){
  return resname;
}

std::string& Atom::getChainId(){
  return chainid;
}

int& Atom::getResId(){
  return resid;
}

std::string& Atom::getICode(){
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

int& Atom::getSel(){
  return sel;
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

void Atom::setAtmName(){
  this->atmname.clear();
}

void Atom::setAlt(const std::string& altin){
  this->alt=altin;
}

void Atom::setAlt(){
  this->alt.clear();
}

void Atom::setResName(const std::string& resnamein){
  this->resname=resnamein;
}

void Atom::setResName(){
  this->resname.clear();
}

void Atom::setChainId(const std::string& chainidin){
  this->chainid=chainidin;
}

void Atom::setChainId(){
  this->chainid.clear();
}

void Atom::setResId(const int& residin){
  this->resid=residin;
}

void Atom::setICode(const std::string& icodein){
  this->icode=icodein;
}

void Atom::setICode(){
  this->icode.clear();
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

void Atom::setSegId(){
  this->segid.clear();
}

void Atom::setSel(const int selin){
  this->sel=selin;
}
