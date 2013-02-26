//Sean M. Law

#ifndef TRAJECTORY_H
#define TRAJECTORY_H

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <climits>

class Trajectory {
  private:
    std::string format;
		char hdr[4];
		int nframe; //ICNTRL[0], NFILE
		double tstart; //ICNTRL[1], NPRIV
		unsigned int natom;
		double tstep;
		double dof;
		bool periodic;
		double pbx;
		double pby;
		double pbz;
		bool fixed;
		unsigned int version;
		std::string endian;
		std::string title;
		std::vector<double> x;
		std::vector<double> y;
		std::vector<double> z;
		typedef union {
			unsigned int ui;
			int i;
			char c[4];
			float f;
			double d;
		} binbuf;

  public:
    bool findFormat(std::ifstream &trajin);
		Trajectory::binbuf* readFortran(std::ifstream &trajin);
		void readHeader(std::ifstream &trajin);
};


#endif
