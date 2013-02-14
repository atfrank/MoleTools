//Sean M. Law

#ifndef MISC_H
#define MISC_H

#include <vector>
#include <string>
#include <sstream>
#include <iostream>
#include <algorithm>

class Misc {
  private:

  public:
    static std::vector<std::string> split (const std::string &str, const std::string &delim, const bool repeat=true);
    static bool isdigit (const std::string &str);
    static bool isalpha (const std::string &str);
    static bool isrange (const std::string &str);
    static std::string trim (const std::string &str, const std::string t=" ");
    static std::string processRange (const std::string &start, const std::string &end);
};

#endif
