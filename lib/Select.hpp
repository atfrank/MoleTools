//Sean M. Law
#ifndef SELECT_H
#define SELECT_H

#include <Molecule.hpp>
#include <map>

class Select {
  private:
    std::map<std::string, int> chainIdMap;
    std::map<std::string, int> segIdMap;
    std::map<std::string, int> resNameMap;
    std::map<std::string, int> resIdMap;
    std::map<std::string, int> atmNameMap;

  public:
    static void makeSel(Molecule* mol, std::string sel);
    void parseSel(std::string sel);
};

#endif
