//Sean M. Law

//Contact Appearance/Disappearance Order (CAO)

#include "Misc.hpp"

#include <iostream>
#include <fstream>
#include <cstdlib>

#define MAXINPFILES 4096
#define MAXLINESIZE 4096

void usage(){
	std::cerr << std::endl << std::endl;
	std::cerr << "Usage:   contactOrder [-options] <file(s)>" << std::endl;
	std::cerr << "Options: [-min contacts]" << std::endl;
	std::cerr << "         [-skip lines]" << std::endl;
  exit(0);
}

int main (int argc, char **argv){

  int i,l,m,n;
	unsigned int j, k, p;
  std::string currArg;
	std::vector<std::string> ifiles;
	std::ifstream inpFile;
	std::istream* inp;
	std::string line;
	int skip;
	std::vector<int> s;
	//Note that contact "order" and contact "rank" mean the same thing 
	//except that order is used for one binding even while rank is used
	//for the normalized sum of all binding events
	std::vector<std::pair<int,int> > order;
	std::vector<std::vector<double> > rank;
	std::pair<std::vector<std::pair<int, int> >::iterator, std::vector<std::pair<int,int> >::iterator> range;
	unsigned int time;
	int nline;
	int nevents;
	bool trackFlag;
	int min; //Not-implemented yet! Adjust the number of contacts needed for "fully" bound state
	int lastTime;
	int start;
	int stop;
	int nrange;

	ifiles.clear();
	skip=0;
	nline=0;
	order.clear();
	time=0;
	nevents=0;
	trackFlag=false;
	min=-1;
	lastTime=-1;

  for (i=1; i<argc; i++){
    currArg=argv[i];
    if (currArg == "-h" || currArg == "-help"){
      usage();
    }
		else if (currArg == "-min"){
			currArg=argv[++i];
			std::stringstream(currArg) >> min;
		}
    else if (currArg == "-skip"){
      currArg=argv[++i];
			std::stringstream(currArg) >> skip;
			if (skip > 0){
				skip--;
			}
    }
    else{
			ifiles.push_back(currArg);
    }
  }

	//Loop through each file
	for (j=0; j< ifiles.size(); j++){
		if (ifiles.at(j) == "-"){
    	inp=&std::cin;
  	}
  	else{
    	inpFile.open((ifiles.at(j)).c_str());
    	inp=&inpFile;
  	}
    std::cerr << "Processing file " << ifiles.at(j) << "..." << std::endl;
		while (inp->good() && !(inp->eof())){
			getline(*inp, line);
			if (line.size() == 0){
				continue;
			}
			if (nline != skip){
				nline++;
				continue;
			}
			else{
				//Assess snapshot
				Misc::splitNum(line, " \t", s, false); //Split on one or more consecutive whitespaces
				if (min >= 0 && min <= s.at(2)){
					s.at(2)=min;
				}
				if (order.size() == 0){
//					order.resize((s.size()-3)/2);
					order.resize(s.at(2));
					for (k=0; k< order.size(); k++){
						order.at(k)=std::make_pair(k,0);
					}
				}
				if (rank.size() == 0){
//					rank.resize((s.size()-3)/2);
					rank.resize(s.at(2));
					for (k=0; k< rank.size(); k++){ //Contact
//						rank.at(k).resize((s.size()-3)/2);
						rank.at(k).resize(s.at(2));
						for (p=0; p< rank.size(); p++){ //Rank
            	rank.at(k).at(p)=0.0;
						}
          }
				}
				if (s.at(0) == 0 && trackFlag == false){
					time=0;
					for (k=0; k< order.size(); k++){
						order.at(k).second=0;
					}	
					trackFlag=true;
				}
				else{
					if (s.at(0) > 0 && trackFlag == true){
						time++;
						l=0;
						for (k=3; k< s.size(); k=k+2){
							if (s.at(k) == 0){
								order.at(l).second=0; //Set time to zero if contact is broken
							}
							else{
								//Contact is formed
								if (order.at(l).second == 0){
									order.at(l).second=time; //Only record time if the contact was broken and now formed
								}
							}
							l++;
						}
						if (s.at(0) == s.at(2)){ //All contacts are formed
            	//Sort by time
							std::sort(order.begin(), order.end(), Misc::sortPairSecond);

							//Find range of contacts for each time and add to rank
							lastTime=-1;
							nrange=1;
							l=0;
							for (k=0; k< order.size(); k=k+nrange){
								if (lastTime != order.at(k).second){
									range=std::equal_range(order.begin()+k, order.end(), order.at(k), Misc::sortPairSecond);

									start = range.first-order.begin(); //Start of last time
									stop = range.second-order.begin(); //Start of next time = Stop of last time +1
									nrange=range.second-range.first; 

									//Add current order to rank
									for (m=start; m< stop; m++){
										for (n=start; n< stop; n++){ //Contact index, actual contact number is "order.at(n).first"
											rank.at(order.at(n).first).at(l)+=1.0/nrange;
										}
										l++;
									}
								}
								lastTime=order.at(k).second;
							}

							//Return order of contacts to original form
							for (k=0; k< rank.size(); k++){
		            order.at(k)=std::make_pair(k,0);
		          }
            	trackFlag=false;
            	nevents++;
          	}
					}
				}
				nline=0;
			}
		}
		if (ifiles.at(j) != "-"){
			inpFile.close();
		}
	}

	//Normalize ranks by nevents
	for (j=0; j< rank.size(); j++){
		for (k=0; k< rank.size(); k++){
			std::cout << j+1 << " " << k+1 << " " << rank.at(j).at(k)/nevents << " " << nevents << std::endl;
		}
		std::cout << std::endl;
	}
  return 0;
}
