//Sean M. Law

#include "Vector.hpp"

#include <iostream>

Vector::Vector(){
  xcoor=0.0;
  ycoor=0.0;
  zcoor=0.0;
}

Vector::Vector(double xcoorin, double ycoorin, double zcoorin){
  xcoor=xcoorin;
  ycoor=ycoorin;
  zcoor=zcoorin;
}

Vector::Vector(const Vector& vec){
  xcoor=vec.xcoor;
  ycoor=vec.ycoor;
  zcoor=vec.zcoor;
}

Vector& Vector::operator= (const Vector& vec){
  xcoor=vec.xcoor;
  ycoor=vec.ycoor;
  zcoor=vec.zcoor;
  return(*this);
}

Vector& Vector::operator= (const double val){
  xcoor=val;
  ycoor=val;
  zcoor=val;
  return(*this);
}

//Addition
Vector& Vector::operator+ (const Vector& vec){
  xcoor+=vec.xcoor;
  ycoor+=vec.ycoor;
  zcoor+=vec.zcoor;
  return(*this);
}

Vector& Vector::operator+= (const Vector& vec){
  xcoor+=vec.xcoor;
  ycoor+=vec.ycoor;
  zcoor+=vec.zcoor;
  return(*this);
}

Vector& Vector::operator+ (const double val){
  xcoor+=val;
  ycoor+=val;
  zcoor+=val;
  return(*this);
}

Vector& Vector::operator+= (const double val){
  xcoor+=val;
  ycoor+=val;
  zcoor+=val;
  return(*this);
}

//Subtraction
Vector& Vector::operator- (const Vector& vec){
  xcoor-=vec.xcoor;
  ycoor-=vec.ycoor;
  zcoor-=vec.zcoor;
  return(*this);
}

Vector& Vector::operator-= (const Vector& vec){
  xcoor-=vec.xcoor;
  ycoor-=vec.ycoor;
  zcoor-=vec.zcoor;
  return(*this);
}

Vector& Vector::operator- (const double val){
  xcoor-=val;
  ycoor-=val;
  zcoor-=val;
  return(*this);
}

Vector& Vector::operator-= (const double val){
  xcoor-=val;
  ycoor-=val;
  zcoor-=val;
  return(*this);
}

//Multiplication
Vector& Vector::operator* (const Vector& vec){
  xcoor*=vec.xcoor;
  ycoor*=vec.ycoor;
  zcoor*=vec.zcoor;
  return(*this);
}

Vector& Vector::operator*= (const Vector& vec){
  xcoor*=vec.xcoor;
  ycoor*=vec.ycoor;
  zcoor*=vec.zcoor;
  return(*this);
}

Vector& Vector::operator* (const double val){
  xcoor*=val;
  ycoor*=val;
  zcoor*=val;
  return(*this);
}

Vector& Vector::operator*= (const double val){
  xcoor*=val;
  ycoor*=val;
  zcoor*=val;
  return(*this);
}

//Division
Vector& Vector::operator/ (const Vector& vec){
  xcoor/=vec.xcoor;
  ycoor/=vec.ycoor;
  zcoor/=vec.zcoor;
  return(*this);
}

Vector& Vector::operator/= (const Vector& vec){
  xcoor/=vec.xcoor;
  ycoor/=vec.ycoor;
  zcoor/=vec.zcoor;
  return(*this);
}

Vector& Vector::operator/ (const double val){
  xcoor/=val;
  ycoor/=val;
  zcoor/=val;
  return(*this);
}

Vector& Vector::operator/= (const double val){
  xcoor/=val;
  ycoor/=val;
  zcoor/=val;
  return(*this);
}

Vector Vector::operator- () const {
  return Vector(-xcoor,-ycoor,-zcoor);
}

Vector Vector::operator- (const Vector& vec) const {
  return Vector(xcoor-vec.xcoor,ycoor-vec.ycoor,zcoor-vec.zcoor);
}

Vector Vector::operator+ (const Vector& vec) const {
  return Vector(xcoor+vec.xcoor,ycoor+vec.ycoor,zcoor+vec.zcoor);
}

double Vector::dot (const Vector& vec) const { //Dot Product
	return xcoor*vec.xcoor+ycoor*vec.ycoor+zcoor*vec.zcoor;
}

Vector Vector::cross (const Vector& vec) const { //Cross Product
  return Vector(ycoor*vec.zcoor - zcoor*vec.ycoor,
                zcoor*vec.xcoor - xcoor*vec.zcoor,
                xcoor*vec.ycoor - ycoor*vec.xcoor);
}

double Vector::norm () const { //Normal
  return sqrt(xcoor*xcoor+ycoor*ycoor+zcoor*zcoor);
}

double Vector::distance (const Vector& u, const Vector& v) {
	double dx, dy, dz;
	dx=u.xcoor-v.xcoor;
	dy=u.ycoor-v.ycoor;
	dz=u.zcoor-v.zcoor;

  return sqrt(dx*dx+dy*dy+dz*dz);
}

double Vector::angle (const Vector& u, const Vector& v, const Vector& w){
  double angle;

  return angle;
}

double Vector::dihedral (const Vector& t, const Vector& u, const Vector& v, const Vector& w) {
  double dihedral;
	Vector dx, dy, dz, p1, p2, p3;
	double np1, np2, dp1, dp2, ts;

	dx=t-u;
	dy=u-v;
	dz=w-v; //This is correct!

	p1=dx.cross(dy);

	np1=p1.norm();
	p1=p1/np1;

	p2=dz.cross(dy);
  np2=p2.norm();
  p2=p2/np2;

	dp1=p1.dot(p2); //Dot product

	ts=1.0-dp1*dp1;
	ts=(ts<0.0)?0.0:sqrt(ts);
	dihedral=PI/2.0-atan2(dp1,ts);

	p3=p1.cross(p2);
	
	dp2=p3.dot(dy); //Dot product

	if (dp2 > 0.0){
		dihedral=-dihedral;
	}

  return dihedral/PI*180.0;
}
