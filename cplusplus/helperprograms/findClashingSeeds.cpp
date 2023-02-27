#include "mstoptions.h"
#include "mstrotlib.h"
#include "mstsystem.h"

#include "generateseeds.h"
#include "residuecontact.h"
#include "residueframe.h"
#include "utilities.h"

int main(int argc, char *argv[]) {
    MstOptions op;
    op.setTitle("Loads a binary file and checks whether the given seed clashes with sidechain/backbone of another structure");
    op.addOption("seedBin","Path to seed binary file. In read mode, will load this file. In write mode will either A) load the file, if one already exists or B) create a new one",true);
    op.addOption("targetPDB","Path to a PDB structure (read-mode only). If provided, will combine each structure from the binary file with this target structure (useful if they form a complex)",false);
    op.setOptions(argc,argv);
        
    string seedBinPath = op.getString("seedBin");
    string targetPDBPath = op.getString("targetPDB","");

    if (!MstSys::fileExists(seedBinPath)) MstUtils::error("When in 'read' mode, must provide the path to an existing seed binary file path");

    seedBinaryFile seedBin(seedBinPath);
    seedBin.scanFilePositions();

    Structure target = Structure(targetPDBPath);
    clashChecker checker(target);

    fstream f_out;
    MstUtils::openFile(f_out,"clashing_seeds.csv",fstream::out);
    f_out << "name,n_clashes" << endl;

    string dirName = "output";
    MstSys::cmkdir(dirName);

    seedBin.reset();
    int n_clashes = 0;
    while (seedBin.hasNext()) {
        Structure* seed = seedBin.next();
        n_clashes = checker.countClashesToStructure(seed->getResidues()); 
        if (n_clashes > 0) {
            cout << "Seed " << seed->getName() << " has " << n_clashes << " clash(es)" << endl;
            f_out << seed->getName() << "," << n_clashes << endl;
            // seed->writePDB(dirName+"/"+seed->getName()+".pdb");
        }
        delete seed;
    }

    cout << "Done!" << endl;
    return 0;
}