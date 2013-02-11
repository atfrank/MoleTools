//Sean M. Law

#include <Misc.hpp>

std::vector<std::string> Misc::split (const std::string &str, const std::string &delim, const bool repeat){
  std::vector<std::string> tokens;
  size_t p0=0;
  size_t p1=std::string::npos;

  //"repeat" = true means that a blank string is added when there are
  //back-to-back delimiters. Otherwise, repeat=false ignores back-to-back delimiters.

  p1=str.find_first_of(delim,p0);

  while (p1 != std::string::npos){
    if (p1-p0 > 1){
      tokens.push_back(str.substr(p0,p1-p0));
    }
    else{
      if (repeat){
        tokens.push_back(str.substr(p0,p1-p0));
      }
    }
    p0=p1+1;
    p1=str.find_first_of(delim, p0);
  }
  //After last delimiter
  if (p1-p0 > 1){
   tokens.push_back(str.substr(p0,p1-p0));
  }
  else{
    if (repeat){
      tokens.push_back(str.substr(p0,p1-p0));
    }
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

