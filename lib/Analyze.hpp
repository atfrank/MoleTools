//Sean M. Law

#include "Molecule.hpp"
#include "Vector.hpp"
#include "Constants.hpp"

class Analyze {
	private:

	public:
		static Vector centerOfGeometry(Molecule* mol, bool selFlag=true);
		static double rmsd (Molecule* cmpmol, Molecule* refmol);
		static double distance (const Vector& u, const Vector& v);
		static double angle (const Vector& u, const Vector& v, const Vector& w);
    static double dihedral (const Vector& t, const Vector& u, const Vector& v, const Vector& w);
		static double distance (Molecule* sel1, Molecule* sel2, bool selFlag=true);
		static double angle (Molecule* sel1, Molecule* sel2, Molecule* sel3, bool selFlag=true);
	  static double dihedral (Molecule* sel1, Molecule* sel2, Molecule* sel3, Molecule* sel4,bool selFlag=true);
};
