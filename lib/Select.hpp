//Sean M. Law
#ifndef SELECT_H
#define SELECT_H

#include <Molecule.hpp>

#include <vector>

class Select {
  private:
    struct Selection {
      bool negAll;
      std::vector<std::string> chainids;
      std::vector<std::string> segids;
      std::vector<std::string> resnames;
      std::vector<int> resids;
      std::vector<std::string> atmnames;
      void clear();
    };
    std::vector<Selection> selVec;

  public:
    static void makeSel(Molecule* mol, std::string selin);
    void parseSel(std::string selin);

};

#endif
