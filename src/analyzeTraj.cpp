//Sean M. Law

#include "Trajectory.hpp"
#include "Molecule.hpp"
#include "Misc.hpp"
#include "Analyze.hpp"

#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <vector>

void usage(){
  std::cerr << std::endl << std::endl;
  std::cerr << "Usage:   analyzeTraj [-options] <-pdb PDBFILE> <TRAJfile(s)>" << std::endl;
  std::cerr << "Options: [-sel selection] [-cog selection]" << std::endl;
  std::cerr << "         [-dsel selection selection ] [-dist selection selection]" << std::endl;
	std::cerr << "         [-tsel selection selection selection] [-angle selection selection selection]" << std::endl;
	std::cerr << "         [-qsel selection selection selection selection]" << std::endl;
	std::cerr << "         [-dihedral selection selection selection selection]" << std::endl;
	std::cerr << "         [-fit fitPDB] [-fitsel selection]" << std::endl;
	std::cerr << "         [-rmsd selection]" << std::endl;
	std::cerr << "         [-rmsf selection]" << std::endl;
  std::cerr << "         [-skip frames] [-start frame] [-stop frame]" << std::endl;
	std::cerr << "         [-average selection]" << std::endl;
  std::cerr << "         [-eda selection covarout | -project selection mode1[:mode2[:...:[modeN]]] covarin]" << std::endl;
  std::cerr << "         [-list file]" << std::endl;
  std::cerr << "         [-verbose]" << std::endl;
	exit(0);
}

int main (int argc, char **argv){

  int i;
	unsigned int itrj, ianalysis, iline, j;
  std::vector<std::string> trajs;
  Molecule *mol;
	Molecule *fitmol; //For fitting only
  std::string pdb;
	std::string refpdb;
	std::string fitpdb;
	bool fit=false;
	std::string fitsel=":.";
  std::string currArg;
  std::ifstream trjin;
	Trajectory *ftrjin;
  int skip=0;
  int start=0;
	int stop=-1;
	bool startFlag=false;
	std::string line;
	std::vector<std::string> s;
	Vector xyz;
  std::string flist;
  std::ifstream listFile;
  std::istream* listinp;
  std::vector<unsigned int> frameList;
  unsigned int listinx;
  unsigned int nline;
  std::vector<unsigned int> sframes;
  unsigned int lastFrame;
  unsigned int nframe;
  bool timeseries;
  std::vector<unsigned int> modes;
  bool verbose;

	std::vector<Analyze *> analyses;
	Analyze *anin;

	iline=1;
  pdb.clear();
	fitpdb.clear();
  mol=NULL;
	fitmol=NULL;
	ftrjin=NULL;
  flist.clear();
  listinx=0;
  lastFrame=0;
  nframe=0;
  verbose=false;
	nline=0;
  timeseries=false;
  modes.clear();


  for (i=1; i<argc; i++){
    currArg=argv[i];
    if (currArg.compare("-h") == 0 || currArg.compare("-help") == 0){
      usage();
    }
    else if (currArg.compare("-pdb") == 0){
      currArg=argv[++i];
      pdb=currArg;
    }
    else if (currArg.compare("-sel") == 0 || currArg.compare("-nsel") == 0 || currArg.compare("-cog") == 0){
			anin=new AnalyzeCOG;
      currArg=argv[++i];
			anin->addSel(currArg);
			analyses.push_back(anin);
      timeseries=true;
    }
    else if (currArg.compare("-dsel") == 0 || currArg.compare("-dist") == 0 || currArg.compare("-distance") == 0){
			anin=new AnalyzeDistance;
      currArg=argv[++i];
			anin->addSel(currArg);
      currArg=argv[++i];
			anin->addSel(currArg);
			analyses.push_back(anin);
      timeseries=true;
    }
		else if (currArg.compare("-tsel") == 0 || currArg.compare("-angle") == 0){
			anin=new AnalyzeAngle;
			currArg=argv[++i];
			anin->addSel(currArg);
			currArg=argv[++i];
			anin->addSel(currArg);
			currArg=argv[++i];
			anin->addSel(currArg);
			analyses.push_back(anin);
      timeseries=true;
		}
		else if (currArg.compare("-qsel") == 0 || currArg.compare("-dihedral") == 0){
			anin=new AnalyzeDihedral;
      currArg=argv[++i];
			anin->addSel(currArg);
      currArg=argv[++i];
			anin->addSel(currArg);
      currArg=argv[++i];
			anin->addSel(currArg);
			currArg=argv[++i];
			anin->addSel(currArg);
			analyses.push_back(anin);
      timeseries=true;
		}
		else if (currArg.compare("-fit") == 0){
			fit=true;
			currArg=argv[++i];
			fitpdb=currArg;
		}
		else if (currArg.compare("-fitsel") == 0){
			fit=true;
			currArg=argv[++i];
			fitsel=currArg;
		}
		else if (currArg.compare("-rmsd") == 0 || currArg.compare("-rms") == 0){
			anin=new AnalyzeRMSD;
			currArg=argv[++i];
			anin->addSel(currArg);
      analyses.push_back(anin);
      timeseries=true;
		}
		else if (currArg.compare("-rmsf") == 0){
			anin=new AnalyzeRMSF;
			currArg=argv[++i];
			anin->addSel(currArg);
			analyses.push_back(anin);
		}
		else if (currArg.compare("-average") == 0){
			anin=new AnalyzeAverage;
			currArg=argv[++i];
			anin->addSel(currArg);
			analyses.push_back(anin);
		}
    else if (currArg.compare("-eda") == 0 || currArg.compare("-essential") == 0){
      anin=new AnalyzeEDA;
      currArg=argv[++i];
      anin->addSel(currArg);
      currArg=argv[++i];
      anin->setOutput(currArg);
      analyses.push_back(anin);
    }
    else if (currArg.compare("-project") == 0 || currArg.compare("-projection") == 0){
      anin=new AnalyzeProjection;
      currArg=argv[++i];
      anin->addSel(currArg);
      currArg=argv[++i];
      Misc::splitNum(currArg, ":", modes, false);
      anin->addModes(modes);
      currArg=argv[++i];
      anin->setInput(currArg);
      analyses.push_back(anin);
    //  timeseries=true;
    }
    else if (currArg.compare("-skip") == 0){
      currArg=argv[++i];
      std::stringstream(currArg) >> skip;
    }
    else if (currArg.compare("-start") == 0){
      currArg=argv[++i];
      std::stringstream(currArg) >> start;
      start--;
			startFlag=true;
    }
		else if (currArg.compare("-stop") == 0){
			currArg=argv[++i];
			std::stringstream(currArg) >> stop;
		}
    else if (currArg.compare("-list") == 0){
      currArg=argv[++i];
      flist=currArg;
    }
    else if (currArg.compare("-verbose") == 0 || currArg.compare("-v") == 0){
      verbose=true;
    }
    else{
      trajs.push_back(currArg);
    }
  }

  if (trajs.size() == 0){
    std::cerr << std::endl << "Error: Please provide an input trajectory file" << std::endl << std::endl;
    usage();
  }

	if (pdb.length() == 0){
		std::cerr << std::endl << "Error: Please provide a PDB file via \"-pdb\"" << std::endl;
    usage();
	}
  else {
    mol=Molecule::readPDB(pdb);
		for (j=0; j< analyses.size(); j++){
			analyses.at(j)->preAnalysis(mol); //Make copies of mol from selection
		}
    mol->selAll();
  }

	if (fit == true){
		mol->select(fitsel);
    if (fitpdb.length() > 0){
      if (fitpdb != pdb){
        std::cerr << std::endl << "Warning: fitpdb (\"" << fitpdb;
        std::cerr << "\") is different from pdb (\"" << pdb << "\")!" << std::endl;
      }
      fitmol=Molecule::readPDB(fitpdb);
      fitmol->select(fitsel);
      fitmol=fitmol->clone(true, false); //Clone and delete original
    }
		else{
			fitmol=mol->clone();
		}
		mol->storeSel("fit");
		mol->selAll();
  }

  //Read reference file with list of frames in first column.
  //The frames should be in sequential order and not the frame number
  //for each trajectory.
  //Thus, if you have two trajectories, each with 10, frames then your
  //reference file list should contain numbers between 1-20 (not 1-10)!
  if (flist.length() > 0){
    listFile.open(flist.c_str(), std::ios::in);
    listinp=&listFile;
    while (listinp->good() && !(listinp->eof())){
      getline(*listinp, line);
      nline++;
      Misc::splitNum(line, " \t", sframes, false);
      if (sframes.size() >= 1 && sframes.at(0) > 0){
        if (sframes.at(0) > lastFrame){
          frameList.push_back(sframes.at(0));
          lastFrame=sframes.at(0);
        }
        else{
          std::cerr << std::endl << "Warning: List of frames were not sequential (line ";
          std::cerr << nline << ")!" << std::endl << std::endl;
        }
      }
    }
    if (listFile.is_open()){
      listFile.close();
    }
  }  

	//Process Trajectories

	for (itrj=0; itrj< trajs.size(); itrj++){
		trjin.open(trajs.at(itrj).c_str(), std::ios::binary);
    if (verbose == true){
      std::cerr << "Processing file \"" << trajs.at(itrj) << "\"..." << std::endl;
    }

		if (trjin.is_open()){
			ftrjin=new Trajectory;
      ftrjin->setMolecule(mol);

      if (ftrjin->findFormat(trjin) == true){
				ftrjin->readHeader(trjin);
				if (stop < 0){
					stop=ftrjin->getNFrame();
				}
				if (skip > 0 && startFlag == false){
					start=skip;
				}
        //Loop through desired frames
				for (i=start; i< ftrjin->getNFrame() && i< stop; i=i+1+skip){
					ftrjin->readFrame(trjin, i);
          nframe++;
					//Fit if needed
					if (fit == true){
						ftrjin->getMolecule()->recallSel("fit");
						ftrjin->getMolecule()->lsqfit(fitmol);
						ftrjin->getMolecule()->selAll();
					}
					//Start analyses
          if (flist.length() > 0){
            if (nframe == frameList.at(listinx)){
              if (timeseries == true){
                std::cout << iline << "  ";
                std::cout << ftrjin->getNPriv()*ftrjin->getTStepPS()/ftrjin->getNSavc()+i*ftrjin->getTStepPS();
              }
              for (ianalysis=0; ianalysis< analyses.size(); ianalysis++){
                analyses.at(ianalysis)->runAnalysis();
              }
              if (timeseries == true){
                std::cout << std::endl;
              }
              iline++;
              listinx++;
              if (listinx == frameList.size()){
                break;
              }
            }
          }
          else{
            if (timeseries == true){
					    std::cout << iline << "  ";
					    std::cout << ftrjin->getNPriv()*ftrjin->getTStepPS()/ftrjin->getNSavc()+i*ftrjin->getTStepPS();
            }
					  for (ianalysis=0; ianalysis< analyses.size(); ianalysis++){
						  analyses.at(ianalysis)->runAnalysis();
					  }
            if (timeseries == true){
					    std::cout << std::endl;
            }
					  iline++;
          }
				}
			}
			else{
				std::cerr << "Warning: Skipping unknown trajectory format \"";
				std::cerr << trajs.at(itrj) << "\"" << std::endl;
			}
			if (ftrjin != NULL){
				delete ftrjin;
			}
		}
		trjin.close();
	}

	for (ianalysis=0; ianalysis< analyses.size(); ianalysis++){
    analyses.at(ianalysis)->postAnalysis();
  }

	//delete mol;

  return 0;
}
