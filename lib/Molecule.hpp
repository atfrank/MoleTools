//Sean M. Law

class Molecule {
  private:
    char recname[6]; //Record name: "ATOM  ", "HETATM"
    int  atmnum; //Atom serial number
    char atmname[6]; //Atom name
    char alt; //Alternate location indicator
    char resname[6]; //Residue name
    char chainid; //Chain identifier
    int  resid; //Residue sequence number
    char icode; //Code for insertion of residues
    Vector coor;
    
}
