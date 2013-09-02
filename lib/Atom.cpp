//Sean M. Law

#include "Atom.hpp"

Atom::Atom(){
	pdbid="NoID";
  recname="ATOM";
  atmnum=0;
  atmname="    ";
  alt=" ";
  resname="   ";
  chainid=" ";
  realid=" ";
  resid=0;
  coor=Vector(0.0, 0.0, 0.0);
  occu=0.0;
  bfac=0.0;
  segid="   ";
  sel=true;
  summary="";
	ss="";
}

Atom::Atom(int atmnumin, std::string atmnamein, std::string resnamein, int residin, Vector coorin, std::string segidin){
	pdbid="NoID";
  recname="ATOM";
  atmnum=atmnumin;
  atmname=atmnamein;
  alt=" ";
  resname=resnamein;
  chainid=" ";
  realid=" ";
  resid=residin;
  coor=coorin;
  occu=0.0;
  bfac=0.0;
  segid=segidin;
  sel=true;
  summary="";
	ss="";
}

void Atom::reset(){
	pdbid="NoID";
  recname="ATOM";
  atmnum=0;
  atmname="    ";
  alt=" ";
  resname="   ";
  chainid=" ";
  realid=" ";
  resid=0;
  coor=Vector(0.0, 0.0, 0.0);
  occu=0.0;
  bfac=0.0;
  segid="    ";
  sel=true;
  summary="";
	ss="";
}

void Atom::clone(Atom* atmin){
	pdbid=atmin->getPdbId();
  recname=atmin->getRecName();
  atmnum=atmin->getAtmNum();
  atmname=atmin->getAtmName();
  alt=atmin->getAlt();
  resname=atmin->getResName();
  chainid=atmin->getChainId();
  realid=atmin->getRealId();
  resid=atmin->getResId();
	icode=atmin->getICode();
  coor=atmin->getCoor();
  occu=atmin->getOccu();
  bfac=atmin->getBFac();
  segid=atmin->getSegId();
  sel=atmin->getSel();
  summary=atmin->getSummary();
	ss=atmin->getSS();
}

//Get atom info
std::string& Atom::getPdbId(){
	return pdbid;
}

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

std::string& Atom::getRealId(){
  return realid;
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

bool& Atom::getSel(){
  return sel;
}

std::string& Atom::getSummary(){
  return summary;
}

std::string& Atom::getSS(){
	return ss;
}

//Set atom info
void Atom::setPdbId(const std::string& pdbidin){
	pdbid=pdbidin;
}

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
  atmname="    ";
}

void Atom::setAlt(const std::string& altin){
  alt=altin;
}

void Atom::setAlt(){
  alt=" ";
}

void Atom::setResName(const std::string& resnamein){
  resname=resnamein;
}

void Atom::setResName(){
  resname=" ";
}

void Atom::setChainId(const std::string& chainidin){
  //Modified chain id if duplicated
  chainid=chainidin;
}

void Atom::setRealId(const std::string& realidin){
  //Original chain id from PDB
  realid=realidin;
}

void Atom::setChainId(){
  chainid=" ";
}

void Atom::setRealId(){
  realid=" ";
}

void Atom::setResId(const int& residin){
  resid=residin;
}

void Atom::setICode(const std::string& icodein){
  icode=icodein;
}

void Atom::setICode(){
  icode=" ";
}

void Atom::setX (const double& xin){
  coor.x()=xin;
}

void Atom::setY (const double& yin){
  coor.y()=yin;
}

void Atom::setZ (const double& zin){
  coor.z()=zin;
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
  segid="    ";
}

void Atom::setSel(const bool selin){
  sel=selin;
}

void Atom::setSummary(const std::string& summaryin){
  summary=summaryin;
}

void Atom::setSS(const std::string& ssin){
	ss=ssin;
}
