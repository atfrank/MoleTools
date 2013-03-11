//Sean M. Law

#include "Molecule.hpp"

class Analyze {
	private:

	public:
		static Vector centerOfGeometry(Molecule *mol, bool selFlag=true);
		static double rmsd (Molecule *cmpmol, Molecule *refmol);
};
