//Sean M. Law

#include "Trajectory.hpp"
#include "Molecule.hpp"
#include "Misc.hpp"

#include <fstream>

void usage(){
  std::cerr << std::endl << std::endl;
  std::cerr << "Usage:   traj2PDB [-options] <-pdb PDBfile> <TRAJfile(s)>" << std::endl;
  std::cerr << "Options: [-sel selection]" << std::endl;
  std::cerr << "         [-out tag] [-outsel selection]" << std::endl;
  std::cerr << "         [-skip frames] [-start frame] [-stop frame]" << std::endl;
	std::cerr << "         [-show]" << std::endl;
  std::cerr << "         [-format type] [-chains]" << std::endl;
  exit(0);
}

int main (int argc, char **argv){

  int i;
	unsigned int j;
  std::vector<std::string> trajs;
	std::stringstream fout;
	std::string outname;
	std::string tag="file";
  Molecule *mol;
  Molecule *outmol;
  std::string pdb;
  std::string sel=""; //To match NATOM in trajectory
  std::string currArg;
  std::string outsel="";
  std::ifstream trjin;
	std::ofstream pdbout;
	Trajectory *ftrjin;
  bool show=false;
	int skip=0;
	int start=0;
  int stop=std::numeric_limits<int>::max();
  unsigned int npdb;
  std::string format;
  bool chnFlag;

  pdb.clear();
  mol=NULL;
  outmol=NULL;
	ftrjin=NULL;
  npdb=1;
	format.clear();

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
    else if (currArg.compare("-out") == 0){
      currArg=argv[++i];
      tag=currArg;
    }
    else if (currArg.compare("-outsel") == 0){
      currArg=argv[++i];
      outsel=currArg;
    }
    else if (currArg.compare("-skip") == 0){
      currArg=argv[++i];
      std::stringstream(currArg) >> skip;
    }
    else if (currArg.compare("-start") == 0){
      currArg=argv[++i];
      std::stringstream(currArg) >> start;
      start--;
    }
    else if (currArg.compare("-stop") == 0){
      currArg=argv[++i];
      std::stringstream(currArg) >> stop;
    }
    else if (currArg.compare("-show") == 0){
      show=true;
    }
    else if (currArg.compare("-format") == 0){
      currArg=argv[++i];
      Misc::toupper(currArg);
      format=currArg;
    }
    else if (currArg.compare("-chains") == 0){
      chnFlag=true;
    }
		else if (currArg.compare(0,1,"-") == 0){
      std::cerr << "Warning: Skipping unknown option \"" << currArg << "\"" << std::endl;
    }
    else{
      trajs.push_back(currArg);
    }
  }

  if (trajs.size() == 0){
    std::cerr << std::endl << "Error: Please provide an input trajectory file" << std::endl << std::endl;
    usage();
  }
  else if (pdb.length() == 0){
    std::cerr << std::endl << "Error: Please provide a PDB file via \"-pdb\"" << std::endl;
    usage();
  }
  else {
    mol=Molecule::readPDB(pdb, format);
    if (sel.length() > 0){
      mol->select(sel); //Selection should match trajectory
      mol=mol->clone(true, false); //Delete original molecule after cloning
    }

    if (outsel.length() > 0){
      mol->select(outsel);
      outmol=mol->copy();
    }
  }

	mol->selAll();

	for (j=0; j< trajs.size(); j++){
		trjin.open(trajs.at(j).c_str(), std::ios::binary);

		if (trjin.is_open()){
			ftrjin=new Trajectory;
      ftrjin->setShow(show);
      ftrjin->setMolecule(mol);

      if (ftrjin->findFormat(trjin) == true){
				ftrjin->readHeader(trjin);
        //Loop through desired frames
				for (i=start; i< ftrjin->getNFrame() && i< stop; i=i+1+skip){
					ftrjin->readFrame(trjin, i);
					fout << tag << "." << npdb << ".pdb";
					outname=fout.str();
					pdbout.open(outname.c_str(), std::ios::out);
          std::cerr << "Writing \"" << outname << "\"..." << std::endl;
					if (outsel.length() > 0){
        		pdbout << outmol->writePDB(true, false);
      		}
      		else{
        		pdbout << ftrjin->getMolecule()->writePDB(true, false, chnFlag);
      		}
					pdbout.close();
					fout.str(std::string()); //Clear fout
          npdb++;
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

  return 0;
}
