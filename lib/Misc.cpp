//Sean M. Law

#include "Misc.hpp"

std::string Misc::single (const std::string &resnamein){
	if (resnamein == "T"){

	}
	
	return "";
}

std::vector<std::string> Misc::split (const std::string &str, const std::string &delim, const bool repeat){
  std::vector<std::string> tokens;
  size_t p0=0;
  size_t p1=std::string::npos;

  //"repeat" = true means that a blank string is added when there are
  //back-to-back delimiters. Otherwise, repeat=false ignores back-to-back delimiters.

  p1=str.find_first_of(delim,p0);

  while (p1 != std::string::npos){
    if (p1-p0 > 0){
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
  if (p1-p0 > 0 && p1 != std::string::npos){
    tokens.push_back(str.substr(p0,p1-p0));
  }
  else{
    if (repeat){
      tokens.push_back(str.substr(p0,p1-p0));
    }
  }

  return tokens;
}

bool Misc::isdigit (const std::string &str){
  return str.find_first_not_of("0123456789") == std::string::npos;
}

bool Misc::isalpha (const std::string &str){
  return str.find_first_not_of("ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz") == std::string::npos;
}

bool Misc::isrange (const std::string &str){

  if (str.find_first_not_of("0123456789-") != std::string::npos){
    return false;
  }
  if (str.find("-") == std::string::npos){
    return false;
  }
  if (str.find_first_not_of("-") != 0 || str.find_last_not_of("-") < str.find("-")){
    return false;
  }

  return true;
}

std::string Misc::trim (const std::string &str, const std::string t){
  std::string out=str;

  size_t pos = out.find_last_not_of(t);
  if (pos != std::string::npos){
    if (out.length() != pos+1){
      out.erase(pos+1);
    }
    pos=out.find_first_not_of(t);
    if (pos != 0){
      out.erase(0, pos);
    }
  }
  else{
    out="";
  }

  return out;
}

std::string Misc::processRange (const std::string &start, const std::string &end){
  std::stringstream ss;
  int i, istart, iend;
  
  std::stringstream(start) >> istart;
  std::stringstream(end) >> iend;

  ss << istart;
  istart++;

  for (i=istart; i<=iend; i++){
    ss << "+" << i;
  }

  return ss.str();
}

void Misc::toupper (std::string &str){
  std::transform(str.begin(), str.end(), str.begin(), ::toupper);
}
