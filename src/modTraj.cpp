//Sean M. Law

#include "Trajectory.hpp"
#include "Molecule.hpp"

#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>

void usage(){
  std::cerr << std::endl << std::endl;
  std::cerr << "Usage:   modTraj [-options] <TRAJfile(s)>" << std::endl;
  std::cerr << "Options: [-pdb PDBfile] [-sel selection]" << std::endl;
  std::cerr << "         [-fit fitPDB] [-fitsel selection]" << std::endl;
	std::cerr << "         [-recenter] [-recsel selection]" << std::endl;
  std::cerr << "         [-out TRAJname] [-outsel selection]" << std::endl;
  std::cerr << "         [-skip frames] [-start frame] [-stop frame]" << std::endl;
	std::cerr << "         [-center]" << std::endl;
	std::cerr << "         [-translate dx dy dz]" << std::endl;
	std::cerr << "         [-rotate r1c1 r1c2 r1c3 r2c1 r2c2 r2c3 r3c1 r3c2 r3c3]" << std::endl;
  std::cerr << "         [-list file]" << std::endl;
	std::cerr << "         [-show | -scan]" << std::endl;
	std::cerr << "         [-setnframe int] [-setnpriv int]" << std::endl;
	std::cerr << "         [-setnsavc int] [-setnstep int]" << std::endl;
	std::cerr << "         [-settstep picoseconds]"  << std::endl;
  std::cerr << "         [-verbose]" << std::endl;
  exit(0);
}

int main (int argc, char **argv){

  int i;
	unsigned int j;
  std::vector<std::string> trajs;
  std::string fout;
	bool out=false;
  Molecule *mol;
  Molecule *outmol;
  Molecule *fitmol;
  Molecule *recmol;
  std::string pdb;
	std::string fitpdb;
  std::string sel=""; //To match NATOM in trajectory
  std::string currArg;
  bool fit=false;
  std::string fitsel=":.";
  bool recenter=false;
  std::string recsel=":.";
	bool translate=false;
	double dx, dy, dz;
	Vector dxyz;
	double r1c1, r1c2, r1c3, r2c1, r2c2, r2c3, r3c1, r3c2, r3c3;
  bool rotate=false;
	bool center=false;
  std::string outsel=":.";
  std::ifstream trjin;
  std::ofstream trjout;
	Trajectory *ftrjin;
	Trajectory *ftrjout;
  bool show=false;
	bool scan=false;
  int skip=0;
  int start=0;
  int stop=std::numeric_limits<int>::max();
  bool startFlag=false; //Tracks if user set start value
  std::string flist;
  std::ifstream listFile;
  std::istream* listinp;
  std::vector<unsigned int> frameList;
  unsigned int listinx;
  std::string line;
  unsigned int nline;
  std::vector<unsigned int> s;
  unsigned int nframe;
  unsigned int lastFrame;
	int setnframe; //ICNTRL[1], Number of frames
	int setnpriv; //ICNTRL[2], Number of previous integration steps
	int setnsavc; //ICNTRL[3], Frequency for saving of frames
	int setnstep; //ICNTRL[4],Number of steps in the run that created this file
	double settstep;
  bool verbose;
	unsigned int currFrames;

  pdb.clear();
	fitpdb.clear();
  mol=NULL;
  outmol=NULL;
  fitmol=NULL;
  recmol=NULL;
	ftrjin=NULL;
	ftrjout=NULL;
	dx=0;
	dy=0;
	dz=0;
  r1c1=0.0;
  r1c2=0.0;
  r1c3=0.0;
  r2c1=0.0;
  r2c2=0.0;
  r2c3=0.0;
  r3c1=0.0;
  r3c2=0.0;
  r3c3=0.0;
  flist.clear();
  listinx=0;
  nframe=0; //This is correct!
	currFrames=0;
  lastFrame=0;
  nline=0;
	setnframe=0;
	setnpriv=0;
	setnsavc=0;
	setnstep=0;
	settstep=0.0;
  verbose=false;

  for (i=1; i<argc; i++){
    currArg=argv[i];
    if (currArg.compare("-h") == 0 || currArg.compare("-help") == 0){
      usage();
    }
    else if (currArg.compare("-pdb") == 0){
      currArg=argv[++i];
      pdb=currArg;
    }
    else if (currArg.compare("-sel") == 0){
      currArg=argv[++i];
      sel=currArg;
    }
    else if (currArg.compare("-fit") == 0){
			currArg=argv[++i];
			fitpdb=currArg;
      fit=true;
    }
    else if (currArg.compare("-fitsel") == 0){
      currArg=argv[++i];
      fitsel=currArg;
      fit=true;
    }
    else if (currArg.compare("-recenter") == 0){
      recenter=true;
    }
    else if (currArg.compare("recsel") == 0){
      currArg=argv[++i];
      recsel=currArg;
      recenter=true;
    }
		else if (currArg.compare("-center") == 0){
			center=true;
		}	
		else if (currArg.compare("-translate") == 0){
      currArg=argv[++i];
      std::stringstream(currArg) >> dx;
      currArg=argv[++i];
      std::stringstream(currArg) >> dy;
      currArg=argv[++i];
      std::stringstream(currArg) >> dz;
      dxyz=Vector(dx, dy, dz);
      translate=true;
    }
    else if (currArg.compare("-rotate") == 0){
      currArg=argv[++i];
      std::stringstream(currArg) >> r1c1;
      currArg=argv[++i];
      std::stringstream(currArg) >> r1c2;
      currArg=argv[++i];
      std::stringstream(currArg) >> r1c3;
      currArg=argv[++i];
      std::stringstream(currArg) >> r2c1;
      currArg=argv[++i];
      std::stringstream(currArg) >> r2c2;
      currArg=argv[++i];
      std::stringstream(currArg) >> r2c3;
      currArg=argv[++i];
      std::stringstream(currArg) >> r3c1;
      currArg=argv[++i];
      std::stringstream(currArg) >> r3c2;
      currArg=argv[++i];
      std::stringstream(currArg) >> r3c3;
      rotate=true;
    }
    else if (currArg.compare("-out") == 0){
      currArg=argv[++i];
      fout=currArg;
			out=true;
    }
    else if (currArg.compare("-outsel") == 0){
      currArg=argv[++i];
      outsel=currArg;
			out=true;
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
    else if (currArg.compare("-show") == 0){
      show=true;
			scan=false;
    }
		else if (currArg.compare("-scan") == 0){
			scan=true;
			show=false;
		}
		else if (currArg.compare("-setnframe") == 0){
			currArg=argv[++i];
		  std::stringstream(currArg) >> setnframe;
		}
		else if (currArg.compare("-setnpriv") == 0){
			currArg=argv[++i];
			std::stringstream(currArg) >> setnpriv;
		}
		else if (currArg.compare("-setnsavc") == 0){
			currArg=argv[++i];
			std::stringstream(currArg) >> setnsavc;
		}
		else if (currArg.compare("-setnstep") == 0){
			currArg=argv[++i];
			std::stringstream(currArg) >> setnstep;
		}
		else if (currArg.compare("-settstep") == 0){
			currArg=argv[++i];
	    std::stringstream(currArg) >> settstep;
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
  else if (out == true && fout.length() == 0){
    std::cerr << std::endl << "Error: Please specify an output trajectory via \"-out\" ";
    std::cerr << std::endl << "when using option \"-outsel\"" << std::endl;
  }

  if (fit == true && recenter == true){
    std::cerr << std::endl << "Error: Options \"-fit\" and \"-recenter\" cannot be used simultaneously";
    std::cerr << std::endl;
    usage();
  }

  if (((out == true && outsel.length() > 0)) && pdb.length() == 0){
    std::cerr << std::endl << "Error: Please provide a PDB file via \"-pdb\"" << std::endl;
    usage();
  }
	else if (pdb.length() == 0 && (recenter == true || fit == true)){
		std::cerr << std::endl << "Error: Please provide a PDB file via \"-pdb\"" << std::endl;
    usage();
	}
  else if (pdb.length() > 0){
    mol=Molecule::readPDB(pdb);
    if (sel.length() > 0){
      mol->select(sel); //Selection should match trajectory
      mol=mol->clone(true, false); //Delete original molecule after cloning
    }

    if (outsel.length() > 0){
      mol->select(outsel);
      outmol=mol->copy();
    }
    if (fit == true){
			if (fitpdb.length() == 0){
				std::cerr << std::endl << "Error: Please provide a PDB file for fitting via \"-fit\"" << std::endl;
    		usage();
			}
			fitmol=Molecule::readPDB(fitpdb);
			fitmol->select(fitsel);
      fitmol=fitmol->clone(true, false); //Delete original after cloning
    }
    if (recenter == true){
      mol->select(recsel);
      recmol=mol->copy();
    }
  }
  else{
    //Do Nothing
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
      Misc::splitNum(line, " \t", s, false);
      if (s.size() >= 1 && s.at(0) > 0){
        if (s.at(0) > lastFrame){
          frameList.push_back(s.at(0));
          lastFrame=s.at(0);
        }
        else{
          std::cerr << std::endl << "Warning: List of frames were not sequential (line ";
          std::cerr << nline << ")!" << std::endl << std::endl;
        }
      }
    }
  }

  if (out == true){
    trjout.open(fout.c_str(), std::ios::binary);
		ftrjout=new Trajectory;
    if (outmol != NULL){
      ftrjout->setMolecule(outmol);
    }
  }

	for (j=0; j< trajs.size(); j++){
    if (verbose == true){
      std::cerr << "Processing file \"" << trajs.at(j) << "\"..." << std::endl;
    }
		trjin.open(trajs.at(j).c_str(), std::ios::binary);
    
		if (trjin.is_open()){
			ftrjin=new Trajectory;
      ftrjin->setShow(show);
			ftrjin->setScan(scan);

      if (pdb.length() > 0){
        ftrjin->setMolecule(mol);
				if (fit == true){
					ftrjin->getMolecule()->select(fitsel);
				}
      }

      if (ftrjin->findFormat(trjin) == true){
				ftrjin->readHeader(trjin);
				if (j == 0 && out == true && trjout.is_open()){
					ftrjout->cloneHeader(ftrjin);
          if (start > 0){
            ftrjout->setNPriv(ftrjout->getNPriv()+start*ftrjout->getNSavc());
          }
          if (skip > 0){
            if (startFlag == false){
              start=skip;
              ftrjout->setNPriv(ftrjout->getNPriv()+start*ftrjout->getNSavc());
            }
            ftrjout->setNSavc(ftrjout->getNSavc()*(skip+1));
          }
          ftrjout->setNFrame(0); //Set nframe depending on number of frames written
          ftrjout->writeHeader(trjout);
				}
        //Loop through desired frames
				for (i=start; i< ftrjin->getNFrame() && i< stop; i=i+1+skip){
					ftrjin->readFrame(trjin, i);
          currFrames++;
          if (pdb.length() >0){
						//Transform molecule only if PDB is provided
            if (fit == true){
							ftrjin->getMolecule()->lsqfit(fitmol);
					  }
					  else if (recenter == true){
							ftrjin->getMolecule()->recenter(recmol);
            }
						else if (center == true){
							ftrjin->getMolecule()->center();
						}
						else if (translate == true){
							ftrjin->getMolecule()->translate(dxyz);
						}
						else if (rotate == true){
							ftrjin->getMolecule()->rotate(r1c1, r1c2, r1c3, r2c1, r2c2, r2c3, r3c1, r3c2, r3c3);
						}
            else{
              //Do nothing
            }
					}
          if (out == true && trjout.is_open()){
            if (flist.length() > 0){
              if (currFrames == frameList.at(listinx)){
                ftrjout->writeFrame(trjout, ftrjin);
                listinx++;
                if (listinx == frameList.size()){
                  break;
                }
              }
            }
            else{
              ftrjout->writeFrame(trjout, ftrjin);
							nframe++;
            }
          }
				}
			}
			else{
				std::cerr << "Warning: Skipping unknown trajectory format \"";
				std::cerr << trajs.at(j) << "\"" << std::endl;
			}
			if (ftrjin != NULL){
				delete ftrjin;
			}
		}
	
		trjin.close();
	}

  if (out == true && trjout.is_open()){
    //Re-write header in case of any changes
		if (setnframe > 0){
			ftrjout->setNFrame(setnframe);
		}
		else{
			ftrjout->setNFrame(nframe);
		}
		if (setnpriv > 0){
			ftrjout->setNPriv(setnpriv);
		}
		if (setnsavc > 0){
			ftrjout->setNSavc(setnsavc);
		}
		if (setnstep > 0){
			ftrjout->setNStep(setnstep);
		}
		if (settstep > 0.0){
			ftrjout->setTStep(settstep);
		}
    ftrjout->writeHeader(trjout);
    trjout.close();
		if (ftrjout != NULL){
			delete ftrjout;
		}
  }

  return 0;
}
