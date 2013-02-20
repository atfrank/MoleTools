//Sean M. Law

#ifndef SABA_H
#define SABA_H

#include <Molecule.hpp>

class SABA {
  private:
    std::vector<std::string> ss;
		std::vector<double> helixDist; //Alpha/310 - i', i'+3
		std::vector<double> helixDihe; //Alpha/310 - i', i'+1, i'+2, i'+3
		std::vector<double> threeDist1; //310 - i'+1, i'+2
		std::vector<double> threeDist2; //310 - i'+1, i'+3
		std::vector<double> betaDist; //Parallel/Anti - i', j'
		std::vector<double> paraDist; //Parallel - i'-1, j'-1
		std::vector<double> antiDist1; //Anti - i'+1, j'-1
		std::vector<double> antiDist2; //Anti - i+1, j (C-alpha)

  public:
		static Molecule* getPseudoCenter(Molecule *mol);
};

#endif
