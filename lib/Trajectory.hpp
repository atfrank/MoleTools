//Sean M. Law

#ifndef TRAJECTORY_H
#define TRAJECTORY_H

#include <fstream>
#include <string>
#include <vector>
#include <climits>

class Trajectory {
  private:
    std::string format;
		unsigned int natom;
		unsigned int nframe;
		double tstart;
		double tend;
		double tstep;
		double dof;
		bool periodic;
		double pbx;
		double pby;
		double pbc;
		bool fixed;
		unsigned int version;
		std::string endian;
		std::string title;
		std::vector<double> x;
		std::vector<double> y;
		std::vector<double> z;

  public:
    void getFormat(std::ifstream &trj);
		std::string readFortran();

};


#endif
