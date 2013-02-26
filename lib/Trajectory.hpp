//Sean M. Law

#ifndef TRAJECTORY_H
#define TRAJECTORY_H

#include "Constants.hpp"

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <climits>
#include <cstdlib>

class Trajectory {
  private:
    std::string format;
		char hdr[4];
		int nframe; //ICNTRL[1], NFILE
		int tstart; //ICNTRL[2], NPRIV
		int first;
		int delta; //ICNTRL[3], NSAVC
		double deltat;
		int dof;
		float tstep; //ICNTRL[10], DELTA
		bool crystal;
		double pbx;
		double pby;
		double pbz;
		bool fixed;
		unsigned int version;
		int natom;
		std::string endian;
		std::string title;
		std::vector<float> x;
		std::vector<float> y;
		std::vector<float> z;
		std::vector<int> fixinx;

		typedef union {
			unsigned int ui;
			int i;
			char c[4];
			float f;
		} binbuf;

  public:
    bool findFormat(std::ifstream &trajin);
		Trajectory::binbuf* readFortran(std::ifstream &trajin, int &length);
		void readHeader(std::ifstream &trajin);
		void readFrame(std::ifstream &trjin);
		int getNFrame();
};


#endif
