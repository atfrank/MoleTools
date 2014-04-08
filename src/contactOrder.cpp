//Sean M. Law

//Contact Appearance/Disappearance Order (CAO)

#include "Misc.hpp"

#include <fstream>
#include <limits>
#include <cmath>
#include <algorithm>

void usage(){
	std::cerr << std::endl << std::endl;
	std::cerr << "Usage:   contactOrder [-options] <file(s)>" << std::endl;
	std::cerr << "Options: [-loss]" << std::endl;
	std::cerr << "         [-start line] [-stop line] [-skip lines]" << std::endl;
  exit(0);
}

int main (int argc, char **argv){

  int i,l;
	unsigned int j, k, m, n, p;
  std::string currArg;
	std::vector<std::string> ifiles;
	std::ifstream inpFile;
	std::istream* inp;
	std::string line;
	unsigned int start, stop, skip;
	std::vector<int> s;
	//Note that contact "order" and contact "rank" mean the same thing 
	//except that order is used for one binding even while rank is used
	//for the normalized sum of all binding events
	std::vector<std::pair<int,int> > order;
	std::vector<std::vector<double> > rank;
	std::pair<std::vector<std::pair<int, int> >::iterator, std::vector<std::pair<int,int> >::iterator> range;
	std::vector<int> nbreak; //Total number of times a contact breaks before being permanently formed
  std::vector<std::vector<int> > breakVec;
  std::vector<double> avgBreak;
  std::vector<double> stddevBreak;
	unsigned int time;
	unsigned int nline;
	int nevents;
	bool trackFlag;
	int lastTime;
	unsigned int tstart;
	unsigned int tstop;
	unsigned int nrange;
	bool lossFlag; //Disappearance of contacts
  unsigned int avgTime; //Tracks average time needed for binding/unbinding after first contact is made/lost

	ifiles.clear();
	start=0;
	stop=std::numeric_limits<unsigned int>::max();
	skip=0;
	nline=0;
	order.clear();
	nbreak.clear();
	time=0;
	nevents=0;
	trackFlag=false;
	lastTime=-1;
	lossFlag=false;
  avgTime=0; 

  for (i=1; i<argc; i++){
    currArg=argv[i];
    if (currArg.compare("-h") == 0 || currArg.compare("-help") == 0){
      usage();
    }
		else if (currArg.compare("-loss") == 0){
			lossFlag=true;
		}
		else if (currArg.compare("-start") == 0){
			currArg=argv[++i];
			std::stringstream(currArg) >> start;
		}
		else if (currArg.compare("-stop") == 0){
			currArg=argv[++i];
			std::stringstream(currArg) >> stop;
		}
    else if (currArg.compare("-skip") == 0){
      currArg=argv[++i];
			std::stringstream(currArg) >> skip;
    }
		else if (currArg.compare(0,1,"-") == 0){
      std::cerr << "Warning: Skipping unknown option \"" << currArg << "\"" << std::endl;
    }
    else{
			ifiles.push_back(currArg);
    }
  }

	//Loop through each file
	for (j=0; j< ifiles.size(); j++){
		if (ifiles.at(j).compare("-") == 0){
    	inp=&std::cin;
  	}
  	else{
    	inpFile.open((ifiles.at(j)).c_str());
    	inp=&inpFile;
  	}
    std::cerr << "Processing file " << ifiles.at(j) << "..." << std::endl;
		while (inp->good() && !(inp->eof())){
			getline(*inp, line);
			nline++;
			if (line.size() == 0 || nline < start || nline > stop){
				continue;
			}
			else if (skip > 0 && nline % (skip+1) !=0){
				continue;
			}
			else{
				//Split snapshot
				Misc::splitNum(line, " \t", s, false); //Split on one or more consecutive whitespaces

				//Initialize once
				if (order.size() == 0){
					order.resize(s.at(2));
					for (k=0; k< order.size(); k++){
						order.at(k)=std::make_pair(k,0);
					}
				}
				if (rank.size() == 0){
					rank.resize(s.at(2));
					for (k=0; k< rank.size(); k++){ //Contact
						rank.at(k).resize(s.at(2));
						for (p=0; p< rank.size(); p++){ //Rank
            	rank.at(k).at(p)=0.0;
						}
          }
				}
				if (nbreak.size() == 0){
					nbreak.resize(s.at(2));
					for (k=0; k< nbreak.size(); k++){
					  nbreak.at(k)=0;
					}
				}

				//Assess snapshot
				if (lossFlag == false && s.at(0) == 0 && trackFlag == false){
					//Reset appearance of contacts
					time=0;
					for (k=0; k< order.size(); k++){
						order.at(k).second=0;
					}	
					trackFlag=true;
				}
				else if (lossFlag == true && s.at(0) == s.at(2) && trackFlag == false){
					//Reset disappearance
					time=0;
          for (k=0; k< order.size(); k++){
            order.at(k).second=0;
          }
          trackFlag=true;
				}
				else{
					if (s.at(0) >= 0 && trackFlag == true){
						time++;
						l=0;
						for (k=3; k< s.size(); k=k+2){
							if (s.at(k) == 0){
								if (lossFlag == false){
									if (order.at(l).second > 0){
										nbreak.at(l)++;
									}
									order.at(l).second=0; //Set time to zero if contact is broken
								}
								if (lossFlag == true && order.at(l).second == 0){
									order.at(l).second=time; //Only record time if contact breaks when REVERSE
								}
							}
							else{
								//Contact is formed
								if (lossFlag == false && order.at(l).second == 0){
									order.at(l).second=time; //Only record time if the contact was broken and now formed
								}
								if (lossFlag == true && order.at(l).second > 0){
									order.at(l).second=0; //Set time to zero if contact is still formed when REVERSE
								}
							}
							l++;
						}
						if ((lossFlag == false && s.at(0) == s.at(2)) || (lossFlag == true && s.at(0) == 0)){ //All contacts are formed/loss

            	//Sort by time
							std::sort(order.begin(), order.end(), static_cast<bool (*)(const std::pair<int,int> &, const std::pair<int,int> &)>(Misc::sortPairSecond));

              //Get accumulate time for later averaging
              avgTime+=order.back().second-order.front().second;
							
							//Find range of contacts for each time and add to rank
							lastTime=-1;
							nrange=1;
							l=0;
							for (k=0; k< order.size(); k=k+nrange){
								if (lastTime != order.at(k).second){
									range=std::equal_range(order.begin()+k, order.end(), order.at(k), static_cast<bool (*)(const std::pair<int,int> &, const std::pair<int,int> &)>(Misc::sortPairSecond));

									tstart = range.first-order.begin(); //Start of last time
									tstop = range.second-order.begin(); //Start of next time = Stop of last time +1
									nrange=range.second-range.first; 

									//Add current order to rank
									for (m=tstart; m< tstop; m++){
										for (n=tstart; n< tstop; n++){ //Contact index, actual contact number is "order.at(n).first"
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
              breakVec.push_back(nbreak);
              for (k=0; k< nbreak.size(); k++){
                nbreak.at(k)=0;
              } 
            	nevents++;
          	}
					}
				}
			}
		}
		if (ifiles.at(j) != "-"){
			inpFile.close();
		}
	}

  //Get break average and standard deviation
  avgBreak.resize(nbreak.size(),0.0);
  stddevBreak.resize(nbreak.size(),0.0);
  for (j=0; j< breakVec.size(); j++){
    for (k=0; k< breakVec.at(j).size(); k++){
      avgBreak.at(k)+=breakVec.at(j).at(k);
    }
  }
  for (j=0; j< avgBreak.size(); j++){
    avgBreak.at(j)/=nevents;
    for (k=0; k< breakVec.size(); k++){
      stddevBreak.at(j)+=pow(breakVec.at(k).at(j)-avgBreak.at(j),2.0);
    }
    if (nevents > 1){
      stddevBreak.at(j)=sqrt(stddevBreak.at(j)/(nevents-1.0));
    }
    else{
      stddevBreak.at(j)=sqrt(stddevBreak.at(j));
    }
  }

	//Normalize ranks by nevents
	for (j=0; j< rank.size(); j++){ //Contact number
		for (k=0; k< rank.size(); k++){ //Rank number (appearance/disappearance order)
			std::cout << j+1 << "   " << k+1 << "   " << rank.at(j).at(k)/nevents << "   " << nevents << "   " << (avgTime*1.0)/nevents << "   " << avgBreak.at(k) << " +/- " << stddevBreak.at(k) << std::endl;
		}
		std::cout << std::endl;
	}
  return 0;
}
