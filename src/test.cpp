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
Compiled Using: g++ larmorca.cpp -o larmorca -Wall -O3 -I$MoleTools/lib -lmoletools -L$MoleTools/lib -I$MoleTools */

#include "Atom.hpp"
#include "Residue.hpp"
#include "Molecule.hpp"
#include "Misc.hpp"
#include "Constants.hpp"
#include "LARMOR.hpp"
#include "Analyze.hpp"
#include "Trajectory.hpp"
#include "Mol2.hpp"


#include <iostream>
#include <ctime>
#include <fstream>
#include <limits>
#include <stdlib.h>
#include <sstream>
#include <iomanip>
#include <math.h>
using namespace std;


int main (int argc, char **argv){

	// write MoleTools program that computes the normal for residue 1 in PDB file: model_mini_20.pdb
  
    // dummy variables (not of significance now)
    vector<string> ligandRingSeles;
    vector<double> ligandRingIntensities;
    ligandRingSeles.clear();
    ligandRingIntensities.clear();
  
    // create an empty molecules and then populate with data from external file
	Molecule *mol=NULL;
	Coor refNormal;
    Coor testCorNorm;
    mol = new Molecule;
	mol->cat(Molecule::readPDB("1R2P_2LPS_1.pdb"));
	
	// create temporary molecules for use later on
	Molecule *tmpmol, *reference, *comparison, *tmpmol2, *tmpmol3;
		
	// instantiate LARMOR
	LARMOR * larm=NULL;
	larm = new LARMOR(ligandRingSeles,ligandRingIntensities,mol);

    // selection of reference base
	mol->select(":1.BASE",false,false);
	reference = mol->copy();

    // molecule selection
	mol->select(":1.~10.0",false,false);
	tmpmol = mol->copy();

    // entire molecule selection
    tmpmol->select(":.BASE", false, false);
	comparison = mol->copy();

    // dummy selections used to compare
    tmpmol2 = comparison->copy();
    tmpmol = tmpmol2->copy();
    tmpmol3 = tmpmol2->copy();

    // Normal reference vector and magnitude
    refNormal = larm->base_plane_normal(reference);
    double refMag = sqrt(refNormal.dot(refNormal));

    // data members and instances
    Residue *res;
    Atom *atm;
    stringstream resid;
    string sele1,sele2;
    string resname, atmname;
    double diffDistance;
    const int DISTANCE_LIMIT = 6;
    const int ANGLE_LIMIT = 6;
    int stackingCount = 0;


	// loop over residues in comparison
    for (unsigned int i = 0; i < comparison->getResVecSize();i++) {
        resid.str(""); //Clear stringstream
        res = comparison->getResidue(i);
        resid << res->getResId();
        resname = res->getResName();
        sele1 = ":"+resid.str()+".BASE";
        cout << sele1<< " " << "...."<< resname << endl;
        tmpmol2->select(sele1, false, false);
        tmpmol3 = tmpmol2->copy();
        testCorNorm = larm->base_plane_normal(tmpmol3);
        cout << " x: "<< testCorNorm.x() << " y: "<< testCorNorm.y() << " z: "<< testCorNorm.z() << endl;

        // get distance between comparison base and ref base
        diffDistance = sqrt((pow(testCorNorm.x(),2) + pow(testCorNorm.y(),2) + pow(testCorNorm.z(),2))
                    - (pow(refNormal.x(),2) + pow(refNormal.y(),2) + pow(refNormal.z(),2)));
        cout << "difference: " << diffDistance;

        // get angle between comparison normal and ref normal
        double tempVectorMag = sqrt(testCorNorm.dot(testCorNorm));
        double tempVectorDotRefVector = testCorNorm.dot(refNormal);
        double angleBtwnTempAndRef = acos(tempVectorDotRefVector/(tempVectorMag * refMag)) * (180.0 / PI);
        cout << endl
             << " temp vector magnitude: " << tempVectorMag
             << endl << " ref magnitude: " << refMag
             << endl << " angle between temp and ref normals: "<< angleBtwnTempAndRef << endl << endl;

        // cout the number of bases within a certain distance and angle from reference
        if(diffDistance < DISTANCE_LIMIT && angleBtwnAndRef < ANGLE_LIMIT){
            stackingCount ++;
        }

    }

    // display number that stacks with reference
    cout << "Stacking Count: " << stackingCount << endl;
 		return 0;
}
