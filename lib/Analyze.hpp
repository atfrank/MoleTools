//Sean M. Law

#ifndef ANALYZE_H
#define ANALYZE_H

#include "Molecule.hpp"
#include "Vector.hpp"
#include "Constants.hpp"
#include "Eigen/Eigenvalues"

//Abstract base class (cannot create instance of it!)
class Analyze {
	private:
    //Since this is an abstract base class
    //all members need to be accessed via a function or passed
    //directly to the analysis function
    std::vector<std::string> sel;
		std::vector<Molecule*> mol;
		std::vector<double> tdata; //Time dependent data, maybe for averaging
    Eigen::MatrixXd avgCovar;
		int ndata; //Total number of datapoints
		bool resel; //Re-do selection for each analysis, not implemented yet

	public:
		Analyze ();
		void addSel(const std::string& selin);
		std::string getSel(const int& element);
		unsigned int getNSel();
    void addMol(Molecule* molin);
		void setMol(const int element, Molecule* molin);
		void resizeNMol(const int sizein);
		Molecule* getMol(const int& element);
		unsigned int getNMol();
		void setNData(const int& ndatain);
		int& getNData();
		std::vector<double>& getTDataVec();
    Eigen::MatrixXd& getAvgCovar();
    void initCovar(const unsigned int& xin, const unsigned int& yin);
		
		//Virtual functions
		virtual void setupMolSel(Molecule* molin);
		virtual void preAnalysis(Molecule* molin); 
		virtual void runAnalysis() =0; //Pure virtual function
		virtual void postAnalysis();

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
    static void averageCovariance (Molecule* cmpmol, Molecule* refmol, Eigen::MatrixXd& covarin, int &ndataIO);
		static void pairwiseDistance(Molecule *mol, std::map<std::pair<Atom*, Atom*>, double>& pdin);
		static void allAnglesDihedrals(Molecule *mol, std::map<Atom*, std::vector<double> >& anglesin);
		static void pcasso(Molecule* mol);
};

//Derived classes

class AnalyzeCOG: public Analyze {
	public:
		void runAnalysis();
};

class AnalyzeRMSD: public Analyze {
	public:
		void preAnalysis(Molecule* molin);
		void runAnalysis();
};

class AnalyzeRMSF: public Analyze {
  public:
		void preAnalysis(Molecule* molin);
    void runAnalysis();
		void postAnalysis();
};

class AnalyzeAverage: public Analyze {
	public:
		void preAnalysis(Molecule* molin);
		void runAnalysis();
		void postAnalysis();
};

class AnalyzeDistance: public Analyze {
  public:
		void setupMolSel(Molecule* molin);
    void runAnalysis();
};

class AnalyzeAngle: public Analyze {
  public:
		void setupMolSel(Molecule* molin);
    void runAnalysis();
};

class AnalyzeDihedral: public Analyze {
	public:
		void setupMolSel(Molecule* molin);
		void runAnalysis();
};

class AnalyzeEDA: public Analyze {
  public:
    void preAnalysis(Molecule* molin);
    void runAnalysis();
    void postAnalysis();
};

#endif
