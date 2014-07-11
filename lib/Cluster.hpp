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

#ifndef CLUSTER_H
#define CLUSTER_H

#include <vector>
#include <map>

class Cluster {
  private:
    unsigned int N;
    std::vector<std::vector <double> > dmatrix; //NxN Distance matrix
    std::vector<unsigned int> cluster; //Store cluster numbers for each N datapoints
    std::map<double, int> inx;

  public:
    Cluster();
    Cluster(unsigned int Nin);
    void readDMatrix(std::string finp);
    void setN(unsigned int Nin);
    unsigned int getN();

    //Virtual Functions
    virtual void preCluster();
    virtual void runCluster()=0;
};

class DPCluster: public Cluster {
  private:
    double dc;
    double dmax;
    double maxDelta;
    std::vector<unsigned int> rho;
    std::vector<double> delta;

  public:
    void preCluster();
    void runCluster();
    void setDc(const double& dcin);
    double getDc();
    void setDmax(const double& dmaxin);
    double getDmax();
    void setMaxDelta(const double& maxdeltain);
    double getMaxDelta();
    void setRho();
    void setDelta();

};

class KCluster: public Cluster {
  private: 

  public:
    void runCluster();
};

#endif
