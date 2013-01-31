//Sean M. Law

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

Vector Atom::getCoor (){
  return coor;
}

void Atom::setCoor (Vector coor){
  this->coor=coor;
}

double Atom::x (){
  return coor.x();
}
