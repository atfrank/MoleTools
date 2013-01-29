//Sean M. Law

#include "Vector.hpp"

Vector::Vector(){
  x=0.0;
  y=0.0;
  z=0.0;
}

Vector::Vector(double xin, double yin, double zin){
  x=xin;
  y=yin;
  z=zin;
}

Vector::Vector(const Vector& vec){
      x=vec.x;
      y=vec.y;
      z=vec.z;
    }

Vector& Vector::operator= (const Vector& vec){
  x=vec.x;
  y=vec.y;
  z=vec.z;
  return(*this);
}

Vector& Vector::operator= (const double val){
  x=val;
  y=val;
  z=val;
  return(*this);
}

//Addition
Vector& Vector::operator+ (const Vector& vec){
  x+=vec.x;
  y+=vec.y;
  z+=vec.z;
  return(*this);
}

Vector& Vector::operator+= (const Vector& vec){
  x+=vec.x;
  y+=vec.y;
  z+=vec.z;
  return(*this);
}

Vector& Vector::operator+ (const double val){
  x+=val;
  y+=val;
  z+=val;
  return(*this);
}

Vector& Vector::operator+= (const double val){
  x+=val;
  y+=val;
  z+=val;
  return(*this);
}

//Subtraction
Vector& Vector::operator- (const Vector& vec){
  x-=vec.x;
  y-=vec.y;
  z-=vec.z;
  return(*this);
}

Vector& Vector::operator-= (const Vector& vec){
  x-=vec.x;
  y-=vec.y;
  z-=vec.z;
  return(*this);
}

Vector& Vector::operator- (const double val){
  x-=val;
  y-=val;
  z-=val;
  return(*this);
}

Vector& Vector::operator-= (const double val){
  x-=val;
  y-=val;
  z-=val;
  return(*this);
}

//Multiplication
Vector& Vector::operator* (const Vector& vec){
  x*=vec.x;
  y*=vec.y;
  z*=vec.z;
  return(*this);
}

Vector& Vector::operator*= (const Vector& vec){
  x*=vec.x;
  y*=vec.y;
  z*=vec.z;
  return(*this);
}

Vector& Vector::operator* (const double val){
  x*=val;
  y*=val;
  z*=val;
  return(*this);
}

Vector& Vector::operator*= (const double val){
  x*=val;
  y*=val;
  z*=val;
  return(*this);
}

//Division
Vector& Vector::operator/ (const Vector& vec){
  x/=vec.x;
  y/=vec.y;
  z/=vec.z;
  return(*this);
}

Vector& Vector::operator/= (const Vector& vec){
  x/=vec.x;
  y/=vec.y;
  z/=vec.z;
  return(*this);
}

Vector& Vector::operator/ (const double val){
  x/=val;
  y/=val;
  z/=val;
  return(*this);
}

Vector& Vector::operator/= (const double val){
  x/=val;
  y/=val;
  z/=val;
  return(*this);
}

Vector Vector::operator- () const {
  return Vector(-x,-y,-z);
}

Vector Vector::operator- (const Vector& vec) const {
    return Vector(x-vec.x,y-vec.y,z-vec.z);
}

Vector Vector::operator+ (const Vector& vec) const {
    return Vector(x+vec.x,y+vec.y,z+vec.z);
}

double Vector::operator* (const Vector& vec) const { //Dot Product
    return x*vec.x+y*vec.y+z*vec.z;
}

Vector Vector::cross (const Vector& vec) const { //Cross Product
  return Vector(y*vec.z - z*vec.y,
                z*vec.x - x*vec.z,
                x*vec.y - y*vec.x);
}

double Vector::norm () const { //Normal
  return sqrt(x*x+y*y+z*z);
}
