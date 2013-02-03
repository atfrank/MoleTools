//Sean M. Law

#include <iomanip>

#include "Atom.hpp"

Atom::Atom(){
  recname="ATOM";
  atmnum=0;
  atmname="    ";
  alt=" ";
  resname="   ";
  chainid=" ";
  resid=0;
  coor=Vector(0.0, 0.0, 0.0);
  occu=0.0;
  bfac=0.0;
  segid="   ";
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
  recname="ATOM";
  atmnum=0;
  atmname="    ";
  alt=" ";
  resname="   ";
  chainid=" ";
  resid=0;
  coor=Vector(0.0, 0.0, 0.0);
  occu=0.0;
  bfac=0.0;
  segid="    ";
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
  recname=recnamein;
}

void Atom::setAtmNum(const int& atmnumin){
  atmnum=atmnumin;
}

void Atom::setAtmName(const std::string& atmnamein){
  atmname=atmnamein;
}

void Atom::setAtmName(){
  atmname.clear();
}

void Atom::setAlt(const std::string& altin){
  alt=altin;
}

void Atom::setAlt(){
  alt.clear();
}

void Atom::setResName(const std::string& resnamein){
  resname=resnamein;
}

void Atom::setResName(){
  resname.clear();
}

void Atom::setChainId(const std::string& chainidin){
  chainid=chainidin;
}

void Atom::setChainId(){
  chainid.clear();
}

void Atom::setResId(const int& residin){
  resid=residin;
}

void Atom::setICode(const std::string& icodein){
  icode=icodein;
}

void Atom::setICode(){
  icode.clear();
}

void Atom::setCoor (const Vector& coorin){
  coor=coorin;
}

void Atom::setOccu(const double& occuin){
  occu=occuin;
}

void Atom::setBFac(const double& bfacin){
  bfac=bfacin;
}

void Atom::setSegId(const std::string& segidin){
  segid=segidin;
}

void Atom::setSegId(){
  segid.clear();
}

void Atom::setSel(const int selin){
  sel=selin;
}
