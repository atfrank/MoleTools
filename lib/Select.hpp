//Sean M. Law
#ifndef SELECT_H
#define SELECT_H

#include <Molecule.hpp>

#include <vector>
#include <sstream>
#include <algorithm>
#include <map>

extern std::map<std::string, std::string> selKeys;

class Select {
  private:

  public:
    static void makeSel(Molecule* mol, std::string selin);
    void parseSel(std::string selin);

    //Recursive Descent Parser (RDP)
    static std::vector<Atom*> recursiveDescentParser (const std::string &str, const std::vector<Atom *> &ref, const std::string &group="");
		static std::string getSelValue(const std::string &key);
		static void initKeys();
};

#endif
