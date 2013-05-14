//Sean M. Law

#include "Molecule.hpp"
#include "Vector.hpp"
#include "Constants.hpp"

class Analyze {
	private:
		std::string type;
    std::vector<std::string> sel;
		std::vector<Molecule*> mol;
		std::vector<double> tdata; //Time dependent data, maybe for averaging
		int ndata; //Total number of datapoints
		bool resel; //Re-do selection for each analysis, not implemented yet

	public:
		Analyze ();
		void setType(const std::string& typein);
		std::string getType();
		void addSel(const std::string& selin);
		std::string getSel(const int& element);
		unsigned int getNSel();
		void setupMolSel(Molecule* molin);
    void addMol(Molecule* molin);
		void setMol(const int element, Molecule* molin);
		void resizeNMol(const int sizein);
		Molecule* getMol(const int& element);
		unsigned int getNMol();
		void runAnalysis();
		void postAnalysis();

		//Analysis functions
		static Vector centerOfGeometry(Molecule* mol, bool selFlag=true);
		static double rmsd (Molecule* cmpmol, Molecule* refmol);
		static void rmsf (Molecule* cmpmol, Molecule* refmol, std::vector<double> &tdataIO, int &ndataIO);
		static void averageMol(Molecule* cmpmol, Molecule* refmol, int &ndataIO);
		static double distance (const Vector& u, const Vector& v);
		static double angle (const Vector& u, const Vector& v, const Vector& w);
    static double dihedral (const Vector& t, const Vector& u, const Vector& v, const Vector& w);
		static double distance (Molecule* sel1, Molecule* sel2, bool selFlag=true);
		static double angle (Molecule* sel1, Molecule* sel2, Molecule* sel3, bool selFlag=true);
	  static double dihedral (Molecule* sel1, Molecule* sel2, Molecule* sel3, Molecule* sel4,bool selFlag=true);
};
