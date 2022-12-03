#include "mstsystem.h"
#include "mstoptions.h"

#include "residueframe.h"


int main(int argc, char *argv[]) {
    MstOptions op;
    op.setTitle("Given a protein structure, defines residue frames and writes out PyMol cgo_arrows file for visualization");
    op.addOption("p", "The structure which will be used to compute contact information",true);
    op.setOptions(argc,argv);
    
    augmentedStructure augS(op.getString("p"));
    string structure_name = MstSys::splitPath(augS.getName(),1);

    augS.writeToFile(structure_name);

    return 0;
}