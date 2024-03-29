//Sean M. Law
//Aaron T. Frank
    
/*
This file is part of MoleTools.

MoleTools is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MoleTools is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MoleTools.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "Coor.hpp"

#include <cmath>

Coor::Coor(){
  xcoor=0.0;
  ycoor=0.0;
  zcoor=0.0;
}

Coor::Coor(double xcoorin, double ycoorin, double zcoorin){
  xcoor=xcoorin;
  ycoor=ycoorin;
  zcoor=zcoorin;
}

Coor::Coor(const Coor& vec){
  xcoor=vec.xcoor;
  ycoor=vec.ycoor;
  zcoor=vec.zcoor;
}

Coor& Coor::operator= (const Coor& vec){
  xcoor=vec.xcoor;
  ycoor=vec.ycoor;
  zcoor=vec.zcoor;
  return(*this);
}

Coor& Coor::operator= (const double val){
  xcoor=val;
  ycoor=val;
  zcoor=val;
  return(*this);
}

//Addition
Coor Coor::operator+ (const Coor& vec) const{
  return Coor(xcoor+vec.xcoor,ycoor+vec.ycoor,zcoor+vec.zcoor);
}

Coor& Coor::operator+= (const Coor& vec){
  xcoor+=vec.xcoor;
  ycoor+=vec.ycoor;
  zcoor+=vec.zcoor;
  return(*this);
}

Coor Coor::operator+ (const double val) const{
  return Coor(xcoor+val,ycoor+val,zcoor+val);
}

Coor& Coor::operator+= (const double val){
  xcoor+=val;
  ycoor+=val;
  zcoor+=val;
  return(*this);
}

//Subtraction
Coor Coor::operator- (const Coor& vec) const{
  return Coor(xcoor-vec.xcoor,ycoor-vec.ycoor,zcoor-vec.zcoor);
}

Coor& Coor::operator-= (const Coor& vec){
  xcoor-=vec.xcoor;
  ycoor-=vec.ycoor;
  zcoor-=vec.zcoor;
  return(*this);
}

Coor Coor::operator- (const double val) const{
  return Coor(xcoor-val,ycoor-val,zcoor-val);
}

Coor& Coor::operator-= (const double val){
  xcoor-=val;
  ycoor-=val;
  zcoor-=val;
  return(*this);
}

//Multiplication
Coor Coor::operator* (const Coor& vec) const{
  return Coor(xcoor*vec.xcoor,ycoor*vec.ycoor,zcoor*vec.zcoor);
}

Coor& Coor::operator*= (const Coor& vec){
  xcoor*=vec.xcoor;
  ycoor*=vec.ycoor;
  zcoor*=vec.zcoor;
  return(*this);
}

Coor Coor::operator* (const double val) const{
  return Coor(xcoor*val,ycoor*val,zcoor*val);
}

Coor& Coor::operator*= (const double val){
  xcoor*=val;
  ycoor*=val;
  zcoor*=val;
  return(*this);
}

//Division
Coor Coor::operator/ (const Coor& vec) const{
  return Coor(xcoor/vec.xcoor,ycoor/vec.ycoor,zcoor/vec.zcoor);
}

Coor& Coor::operator/= (const Coor& vec){
  xcoor/=vec.xcoor;
  ycoor/=vec.ycoor;
  zcoor/=vec.zcoor;
  return(*this);
}

Coor Coor::operator/ (const double val) const{
  return Coor(xcoor/val,ycoor/val,zcoor/val);
}

Coor& Coor::operator/= (const double val){
  xcoor/=val;
  ycoor/=val;
  zcoor/=val;
  return(*this);
}


Coor Coor::operator- () const {
  return Coor(-xcoor,-ycoor,-zcoor);
}

double Coor::dot (const Coor& vec) const { //Dot Product
  return xcoor*vec.xcoor+ycoor*vec.ycoor+zcoor*vec.zcoor;
}

Coor Coor::cross (const Coor& vec) const { //Cross Product
  return Coor(ycoor*vec.zcoor - zcoor*vec.ycoor,
                zcoor*vec.xcoor - xcoor*vec.zcoor,
                xcoor*vec.ycoor - ycoor*vec.xcoor);
}

double Coor::norm () const { //Normal
  return sqrt(xcoor*xcoor+ycoor*ycoor+zcoor*zcoor);
}

