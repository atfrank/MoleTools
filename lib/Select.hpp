//Sean M. Law
#ifndef SELECT_H
#define SELECT_H

#include <Molecule.hpp>

#include <vector>

class Select {
  private:
    struct Selection {
      bool negAll;
      std::vector<bool> negChainIds;
      std::vector<std::string> chainids;
      std::vector<bool> negSegIds;
      std::vector<std::string> segids;
      std::vector<bool> negResNames;
      std::vector<std::string> resnames;
      std::vector<bool> negResIds;
      std::vector<int> resids;
      std::vector<bool> negAtmNames;
      std::vector<std::string> atmnames;
      void clear();
    };
    std::vector<Selection> selVec;

  public:
    static void makeSel(Molecule* mol, std::string selin);
    void parseSel(std::string selin);

    //Recursive Descent Parser (RDP)
    std::vector<Atom*> chainRDP (const std::string &str, const std::vector<Atom *> &ref);

    static void makeSelOld(Molecule* mol, std::string selin);
    void parseSelOld(std::string selin);

};

#endif
