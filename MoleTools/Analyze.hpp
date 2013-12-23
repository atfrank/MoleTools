//Sean M. Law

#ifndef ANALYZE_H
#define ANALYZE_H

#include "Enum.hpp"
#include "Molecule.hpp"
#include "Vector.hpp"
#include "Constants.hpp"
#include "Eigen/Eigenvalues"
#include "Eigen/LU"

#include <iostream>
#include <fstream>
#include <algorithm>

//Abstract base class (cannot create instance of it!)
class Analyze {
	private:
    //Since this is an abstract base class
    //all members need to be accessed via a function or passed
    //directly to the analysis function
    std::vector<std::string> sel;
		std::vector<Molecule*> mol;
		std::vector<double> tdata; //Time dependent data, maybe for averaging
		std::vector<std::vector<double> > fdata; //Frame data, cleared after each frame
    Eigen::MatrixXd avgCovar;
    std::vector<unsigned int> modes;
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigen;
		int ndata; //Total number of datapoints
		bool resel; //Re-do selection for each analysis, not implemented yet
    std::string ifile;
    std::string ofile;

	public:
		Analyze ();
		virtual ~Analyze ();
		void addSel(const std::string& selin);
		std::string getSel(const int& element);
		unsigned int getNSel();
    void addMol(Molecule* molin);
		void setMol(const int& element, Molecule* molin);
		void resizeNMol(const int sizein);
		Molecule* getMol(const int& element);
		unsigned int getNMol();
		void setNData(const int& ndatain);
		int& getNData();
		std::vector<double>& getTDataVec();
		std::vector<std::vector<double> >& getFDataVec();
    void setInput(const std::string& fin);
    std::string getInput();
    void setOutput (const std::string& fin);
    std::string getOutput();
    Eigen::MatrixXd& getCovar();
    void addModes(const std::vector<unsigned int>& modesin);
    std::vector<unsigned int>& getModes();
    void initCovar(const unsigned int& xin, const unsigned int& yin);
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> diagonalizeCovar();
    void setEigen(const Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd>& eigenin);
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd>& getEigen();
    void setEigenMode(const unsigned int& modein);
    void readCovariance();
    void writeEigenOverlap(Analyze* cmpin, std::vector<unsigned int>& modein);
		
		//Virtual functions
    virtual void readTopology(Molecule* molin, std::string topin="");
		virtual void setupMolSel(Molecule* molin);
		virtual void preAnalysis(Molecule* molin, std::string topin=""); 
    virtual void preAnalysis();
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
    static std::vector<double> projectModes(Molecule* cmpmol, Molecule* refmol, const Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd>& eigenin, const std::vector<unsigned int>& modesin);
		static void pairwiseDistance(Molecule *mol, std::map<std::pair<Atom*, Atom*>, double>& pdin);
		static void allAnglesDihedrals(Molecule *mol, std::map<Atom*, std::vector<double> >& anglesin);
		static void pcasso(Molecule* mol, std::vector<std::vector<double> > &fdataIO);
    //static std::vector<double> gyration(Molecule* mol);
    Eigen::Matrix3d gyrationTensor(Molecule* mol);
    static double quasiharmonicEntropy(Molecule* mol, const Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd>& eigenin, const std::vector<unsigned int> modesin, double temp=300);
    static double configurationalEntropy(const Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd>& eigenin, const std::vector<unsigned int>& modesin, double cutoffin=1E-10);
};

//Derived classes

class AnalyzeCOG: public Analyze {
	public:
		void runAnalysis();
};

class AnalyzeRMSD: public Analyze {
	public:
		void preAnalysis(Molecule* molin, std::string topin="");
		void runAnalysis();
};

class AnalyzeRMSF: public Analyze {
  public:
		void preAnalysis(Molecule* molin, std::string topin="");
    void runAnalysis();
		void postAnalysis();
};

class AnalyzeAverage: public Analyze {
	public:
		void preAnalysis(Molecule* molin, std::string topin="");
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

class AnalyzeCovariance: public Analyze {
  public:
    void preAnalysis(Molecule* molin, std::string topin="");
    void preAnalysis();
    void runAnalysis();
    void postAnalysis();
};

class AnalyzeProjection: public Analyze {
  public:
    void preAnalysis(Molecule* molin, std::string topin="");
    void runAnalysis();
};

class AnalyzeGyrationTensor: public Analyze {
  public:
    void runAnalysis();
};

class AnalyzeRadiusOfGyration: public Analyze {
  public:
    void runAnalysis();
};

class AnalyzeEllipsoid: public Analyze {
  public:
    void runAnalysis();
};

class AnalyzePairwiseDistance: public Analyze {
  public:
    void runAnalysis();
};

class AnalyzePcasso: public Analyze {
	public:
		void preAnalysis(Molecule* molin, std::string fin="");
		void runAnalysis();
		void printPcasso(PcassoOutEnum out=FEATURES);
};

#endif
