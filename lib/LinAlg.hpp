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

#ifndef LINALG_H
#define LINALG_H

#include <vector>

class LinAlg {
  private:
    std::vector< std::vector<double> > U;
    std::vector< std::vector<double> > V;
    std::vector<double> s;
    
  public:
    void SVD (const std::vector< std::vector <double> > &M);
    void transpose ();
};

#endif
