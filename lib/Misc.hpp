//Sean M. Law

#ifndef MISC_H
#define MISC_H

#include <vector>
#include <string>
#include <sstream>
#include <iostream>

class Misc {
  private:

  public:
    static std::vector<std::string> split (const std::string &str, const std::string &delim, const bool repeat=true);
    static std::vector<int> splitAtoI (const std::string &str, const std::string &delim);
    static std::vector<double> splitAtoF (const std::string &str, const std::string &delim);
};

#endif
