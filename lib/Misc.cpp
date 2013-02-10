//Sean M. Law

#include <Misc.hpp>

std::vector<std::string> Misc::split (const std::string &str, const std::string &delim){
  std::vector<std::string> tokens;
  size_t p0=0;
  size_t p1=std::string::npos;

  while(p0 != std::string::npos){
    p1=str.find_first_of(delim, p0);
    if (p1 != 0){
      std::string token=str.substr(p0, p1-p0);
      tokens.push_back(token);
    }
    else{
      //Nothing inbetween delimiters
      tokens.push_back("");
    }
    p0=str.find_first_not_of(delim, p1);
  }

  return tokens;
}

std::vector<int> Misc::splitAtoI (const std::string &str, const std::string &delim){
  std::vector<int> outTokens;
  std::vector<std::string> inTokens;
  int atoi;

  inTokens=Misc::split(str,delim);
  for (unsigned int i=0; i< inTokens.size(); i++){
    std::stringstream(inTokens.at(i)) >> atoi;
    outTokens.push_back(atoi);
  }

  return outTokens;
}

std::vector<double> Misc::splitAtoF (const std::string &str, const std::string &delim){
  std::vector<double> outTokens;
  std::vector<std::string> inTokens;
  int atof;

  inTokens=Misc::split(str,delim);
  for (unsigned int i=0; i< inTokens.size(); i++){
    std::stringstream(inTokens.at(i)) >> atof;
    outTokens.push_back(atof);
  }

  return outTokens;
}

