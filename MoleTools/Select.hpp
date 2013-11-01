//Sean M. Law
#ifndef SELECT_H
#define SELECT_H

#include "Molecule.hpp"
#include "Misc.hpp"

#include <vector>
#include <sstream>
#include <algorithm>
#include <map>

class Select {
  private:
		std::map<std::string, std::string> selKeysAtm;
		std::map<std::string, std::string> selKeysRes;

  public:
    static void makeSel(Molecule* mol, std::string selin);
    void parseSel(std::string selin);

    //Recursive Descent Parser (RDP)
    std::vector<Atom*> recursiveDescentParser (const std::string &str, const std::vector<Atom *> &ref, const std::string &group="");
		static std::string getSelValue(const std::string &key);
		void initKeys(Molecule *mol);
		bool heavy(const std::string &str, const std::vector<std::string> &heavyVec);
		bool atom(const std::string &str, std::string typein, const std::vector<std::string> &AtomVec);
};

#endif
