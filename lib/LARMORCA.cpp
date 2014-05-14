//Aaron T. Frank
//Sean M. Law

    
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

//Code generated using: awk '{print "{\""$1":"$2":"$3"\","$4"},"}' larmorD_both.dat | tr '\n' ' '

#include "LARMORD.hpp"
#include "LARMORCA.hpp"
#include "Molecule.hpp"
#include "Misc.hpp"

#include <fstream>
#include <stdlib.h>


LARMORCA::LARMORCA (Molecule *mol, const std::string fchemshift){
    /* set random coil chemical shifts */
    this->initializeRandomShifts();

    /* load chemical shifts from file */
    if (fchemshift.length() > 0){
        this->loadCSFile(fchemshift);
    }
}

void LARMORCA::loadCSFile(const std::string fchemshift){
    std::ifstream csFile;
    std::istream* csinp;
    std::string line;
    std::vector<std::string> s;
    std::string n;
    if (fchemshift.length() > 0){
        csFile.open(fchemshift.c_str(), std::ios::in);
        csinp=&csFile;
        while (csinp->good() && !(csinp->eof())){
            getline(*csinp, line);
            Misc::splitStr(line, " ", s, true);
            if (s.size() >= 4){
                n = Misc::trim(s.at(2));
                if( n=="H" || n=="HA" || n=="C" ||n=="CA" || n=="CB" || n=="N" ){
                    this->experimentalCS.insert(std::pair<std::string,double>(Misc::trim(s.at(1))+":"+Misc::trim(s.at(2)),atof(Misc::trim(s.at(3)).c_str())));
                }
            }
        }
    }
}


int LARMORCA::getNShiftAtoms(){
    return (this->shiftAtoms.size());
}

double LARMORCA::getRandomShift(const std::string &key){
    if (this->randomShifts.find (key) == this->randomShifts.end()){
        return 0.0;
    } else {
        return (this->randomShifts.at(key));
    }
}

double LARMORCA::getExperimentalCS(const std::string &key){
    if (this->experimentalCS.find (key) != this->experimentalCS.end()){
        return this->experimentalCS.at(key);
    } else {
        return 0.0;
    }
}

void LARMORCA::initializeRandomShifts(){
    /* generated using: awk '{print "this->randomShifts.insert(std::pair<std::string,double>(\""$1":"$2"\","$3"));"}' randomcoil.dat */
    this->randomShifts.insert(std::pair<std::string,double>("C:ILE",176.4));
    this->randomShifts.insert(std::pair<std::string,double>("C:GLN",176.0));
    this->randomShifts.insert(std::pair<std::string,double>("C:GLY",174.9));
    this->randomShifts.insert(std::pair<std::string,double>("C:GLU",176.6));
    this->randomShifts.insert(std::pair<std::string,double>("C:CYS",174.6));
    this->randomShifts.insert(std::pair<std::string,double>("C:ASP",176.3));
    this->randomShifts.insert(std::pair<std::string,double>("C:SER",174.6));
    this->randomShifts.insert(std::pair<std::string,double>("C:LYS",176.6));
    this->randomShifts.insert(std::pair<std::string,double>("C:PRO",177.3));
    this->randomShifts.insert(std::pair<std::string,double>("C:HID",174.1));
    this->randomShifts.insert(std::pair<std::string,double>("C:HIE",174.1));
    this->randomShifts.insert(std::pair<std::string,double>("C:ASN",175.2));
    this->randomShifts.insert(std::pair<std::string,double>("C:VAL",176.3));
    this->randomShifts.insert(std::pair<std::string,double>("C:THR",174.7));
    this->randomShifts.insert(std::pair<std::string,double>("C:HIS",174.1));
    this->randomShifts.insert(std::pair<std::string,double>("C:TRP",176.1));
    this->randomShifts.insert(std::pair<std::string,double>("C:PHE",175.8));
    this->randomShifts.insert(std::pair<std::string,double>("C:ALA",177.8));
    this->randomShifts.insert(std::pair<std::string,double>("C:MET",173.6));
    this->randomShifts.insert(std::pair<std::string,double>("C:LEU",177.6));
    this->randomShifts.insert(std::pair<std::string,double>("C:ARG",176.3));
    this->randomShifts.insert(std::pair<std::string,double>("C:TYR",175.9));
    this->randomShifts.insert(std::pair<std::string,double>("CB:ILE",38.8));
    this->randomShifts.insert(std::pair<std::string,double>("CB:GLN",29.4));
    this->randomShifts.insert(std::pair<std::string,double>("CB:GLY",999));
    this->randomShifts.insert(std::pair<std::string,double>("CB:GLU",29.9));
    this->randomShifts.insert(std::pair<std::string,double>("CB:CYS",34.5));
    this->randomShifts.insert(std::pair<std::string,double>("CB:ASP",41.1));
    this->randomShifts.insert(std::pair<std::string,double>("CB:SER",63.8));
    this->randomShifts.insert(std::pair<std::string,double>("CB:LYS",33.1));
    this->randomShifts.insert(std::pair<std::string,double>("CB:PRO",33.3));
    this->randomShifts.insert(std::pair<std::string,double>("CB:HID",29.0));
    this->randomShifts.insert(std::pair<std::string,double>("CB:HIE",29.0));
    this->randomShifts.insert(std::pair<std::string,double>("CB:ASN",38.9));
    this->randomShifts.insert(std::pair<std::string,double>("CB:VAL",32.9));
    this->randomShifts.insert(std::pair<std::string,double>("CB:THR",69.8));
    this->randomShifts.insert(std::pair<std::string,double>("CB:HIS",29.0));
    this->randomShifts.insert(std::pair<std::string,double>("CB:TRP",29.6));
    this->randomShifts.insert(std::pair<std::string,double>("CB:PHE",39.6));
    this->randomShifts.insert(std::pair<std::string,double>("CB:ALA",19.1));
    this->randomShifts.insert(std::pair<std::string,double>("CB:MET",32.9));
    this->randomShifts.insert(std::pair<std::string,double>("CB:LEU",42.4));
    this->randomShifts.insert(std::pair<std::string,double>("CB:ARG",30.9));
    this->randomShifts.insert(std::pair<std::string,double>("CB:TYR",37.8));
    this->randomShifts.insert(std::pair<std::string,double>("CA:ILE",61.1));
    this->randomShifts.insert(std::pair<std::string,double>("CA:GLN",55.7));
    this->randomShifts.insert(std::pair<std::string,double>("CA:GLY",45.1));
    this->randomShifts.insert(std::pair<std::string,double>("CA:GLU",56.6));
    this->randomShifts.insert(std::pair<std::string,double>("CA:CYS",56.8));
    this->randomShifts.insert(std::pair<std::string,double>("CA:ASP",54.2));
    this->randomShifts.insert(std::pair<std::string,double>("CA:SER",58.3));
    this->randomShifts.insert(std::pair<std::string,double>("CA:LYS",56.2));
    this->randomShifts.insert(std::pair<std::string,double>("CA:PRO",58.05));
    this->randomShifts.insert(std::pair<std::string,double>("CA:HID",55.0));
    this->randomShifts.insert(std::pair<std::string,double>("CA:HIE",55.0));
    this->randomShifts.insert(std::pair<std::string,double>("CA:ASN",53.1));
    this->randomShifts.insert(std::pair<std::string,double>("CA:VAL",62.2));
    this->randomShifts.insert(std::pair<std::string,double>("CA:THR",61.8));
    this->randomShifts.insert(std::pair<std::string,double>("CA:HIS",55.0));
    this->randomShifts.insert(std::pair<std::string,double>("CA:TRP",57.5));
    this->randomShifts.insert(std::pair<std::string,double>("CA:PHE",57.7));
    this->randomShifts.insert(std::pair<std::string,double>("CA:ALA",52.5));
    this->randomShifts.insert(std::pair<std::string,double>("CA:MET",55.4));
    this->randomShifts.insert(std::pair<std::string,double>("CA:LEU",55.1));
    this->randomShifts.insert(std::pair<std::string,double>("CA:ARG",56.0));
    this->randomShifts.insert(std::pair<std::string,double>("CA:TYR",57.9));
    this->randomShifts.insert(std::pair<std::string,double>("N:ILE",119.9));
    this->randomShifts.insert(std::pair<std::string,double>("N:GLN",119.8));
    this->randomShifts.insert(std::pair<std::string,double>("N:GLY",108.8));
    this->randomShifts.insert(std::pair<std::string,double>("N:GLU",120.2));
    this->randomShifts.insert(std::pair<std::string,double>("N:CYS",118.7));
    this->randomShifts.insert(std::pair<std::string,double>("N:ASP",120.4));
    this->randomShifts.insert(std::pair<std::string,double>("N:SER",115.7));
    this->randomShifts.insert(std::pair<std::string,double>("N:LYS",120.4));
    this->randomShifts.insert(std::pair<std::string,double>("N:PRO",999));
    this->randomShifts.insert(std::pair<std::string,double>("N:HID",118.2));
    this->randomShifts.insert(std::pair<std::string,double>("N:HIE",118.2));
    this->randomShifts.insert(std::pair<std::string,double>("N:ASN",118.7));
    this->randomShifts.insert(std::pair<std::string,double>("N:VAL",119.2));
    this->randomShifts.insert(std::pair<std::string,double>("N:THR",113.6));
    this->randomShifts.insert(std::pair<std::string,double>("N:HIS",118.2));
    this->randomShifts.insert(std::pair<std::string,double>("N:TRP",121.3));
    this->randomShifts.insert(std::pair<std::string,double>("N:PHE",120.3));
    this->randomShifts.insert(std::pair<std::string,double>("N:ALA",123.8));
    this->randomShifts.insert(std::pair<std::string,double>("N:MET",119.6));
    this->randomShifts.insert(std::pair<std::string,double>("N:LEU",121.8));
    this->randomShifts.insert(std::pair<std::string,double>("N:ARG",120.5));
    this->randomShifts.insert(std::pair<std::string,double>("N:TYR",120.3));
    this->randomShifts.insert(std::pair<std::string,double>("H:ILE",8.0));
    this->randomShifts.insert(std::pair<std::string,double>("H:GLN",8.32));
    this->randomShifts.insert(std::pair<std::string,double>("H:GLY",8.3));
    this->randomShifts.insert(std::pair<std::string,double>("H:GLU",8.42));
    this->randomShifts.insert(std::pair<std::string,double>("H:CYS",8.375));
    this->randomShifts.insert(std::pair<std::string,double>("H:ASP",8.34));
    this->randomShifts.insert(std::pair<std::string,double>("H:SER",8.31));
    this->randomShifts.insert(std::pair<std::string,double>("H:LYS",8.29));
    this->randomShifts.insert(std::pair<std::string,double>("H:PRO",999));
    this->randomShifts.insert(std::pair<std::string,double>("H:HID",8.42));
    this->randomShifts.insert(std::pair<std::string,double>("H:HIE",8.42));
    this->randomShifts.insert(std::pair<std::string,double>("H:ASN",8.4));
    this->randomShifts.insert(std::pair<std::string,double>("H:VAL",8.03));
    this->randomShifts.insert(std::pair<std::string,double>("H:THR",8.15));
    this->randomShifts.insert(std::pair<std::string,double>("H:HIS",8.42));
    this->randomShifts.insert(std::pair<std::string,double>("H:TRP",8.25));
    this->randomShifts.insert(std::pair<std::string,double>("H:PHE",8.3));
    this->randomShifts.insert(std::pair<std::string,double>("H:ALA",8.24));
    this->randomShifts.insert(std::pair<std::string,double>("H:MET",8.28));
    this->randomShifts.insert(std::pair<std::string,double>("H:LEU",8.16));
    this->randomShifts.insert(std::pair<std::string,double>("H:ARG",8.23));
    this->randomShifts.insert(std::pair<std::string,double>("H:TYR",8.12));
    this->randomShifts.insert(std::pair<std::string,double>("HA:ILE",4.17));
    this->randomShifts.insert(std::pair<std::string,double>("HA:GLN",4.34));
    this->randomShifts.insert(std::pair<std::string,double>("HA:GLY",3.96));
    this->randomShifts.insert(std::pair<std::string,double>("HA:GLU",4.35));
    this->randomShifts.insert(std::pair<std::string,double>("HA:CYS",4.55));
    this->randomShifts.insert(std::pair<std::string,double>("HA:ASP",4.64));
    this->randomShifts.insert(std::pair<std::string,double>("HA:SER",4.47));
    this->randomShifts.insert(std::pair<std::string,double>("HA:LYS",4.32));
    this->randomShifts.insert(std::pair<std::string,double>("HA:PRO",4.42));
    this->randomShifts.insert(std::pair<std::string,double>("HA:HID",4.73));
    this->randomShifts.insert(std::pair<std::string,double>("HA:HIE",4.73));
    this->randomShifts.insert(std::pair<std::string,double>("HA:ASN",4.74));
    this->randomShifts.insert(std::pair<std::string,double>("HA:VAL",4.12));
    this->randomShifts.insert(std::pair<std::string,double>("HA:THR",4.35));
    this->randomShifts.insert(std::pair<std::string,double>("HA:HIS",4.73));
    this->randomShifts.insert(std::pair<std::string,double>("HA:TRP",4.66));
    this->randomShifts.insert(std::pair<std::string,double>("HA:PHE",4.62));
    this->randomShifts.insert(std::pair<std::string,double>("HA:ALA",4.32));
    this->randomShifts.insert(std::pair<std::string,double>("HA:MET",4.48));
    this->randomShifts.insert(std::pair<std::string,double>("HA:LEU",4.34));
    this->randomShifts.insert(std::pair<std::string,double>("HA:ARG",4.34));
    this->randomShifts.insert(std::pair<std::string,double>("HA:TYR",4.55));
}
