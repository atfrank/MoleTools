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

#ifndef MOL2_H
#define MOL2_H

#include <map>
#include <string>

class Molecule;
class Atom;

class Mol2 {
  private:
    std::map<std::string, int> chnMap;

  public:
    Mol2();
    static Molecule* readMol2 (const std::string ifile, const std::string format="", bool stopFlag=false);
    Atom* processAtomLine (const std::string line, Atom* lastAtom);
};

#endif
