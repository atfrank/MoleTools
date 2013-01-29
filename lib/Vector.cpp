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

//Vector Vector::operator-() const {

//}
