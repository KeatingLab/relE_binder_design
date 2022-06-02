#include "mstoptions.h"
#include "mstrotlib.h"
#include "mstsystem.h"

#include "generateseeds.h"
#include "residuecontact.h"
#include "residueframe.h"
#include "utilities.h"

int main(int argc, char *argv[]) {
    MstOptions op;
    op.setTitle("Loads a binary file, extracts the selected structure, and writes as a PDB file");
    op.addOption("seedBin","Path to seed binary file",true);
    op.addOption("structureName","The name of the structure that is to be extracted from the binary file",false);
    op.addOption("structureList","The path to a file where each line is the name of a structure to be extracted from the binary file",false);
    op.addOption("targetPDB","Path to a PDB structure. If provided, will combine each structure from the binary file with this target structure (useful if they form a complex)",false);
    op.setOptions(argc,argv);

    if ((!op.isGiven("structureName"))&(!op.isGiven("structureList"))) MstUtils::error("Must provide the name of structure, or a list of names");

    string seedBinPath = op.getString("seedBin");
    string structureName = op.getString("structureName","");
    string structureListPath = op.getString("structureList","");
    string targetPDBPath = op.getString("targetPDB","");

    seedBinaryFile seedBin(seedBinPath);
    seedBin.scanFilePositions();

    Structure target;
    vector<Chain*> targetChains;
    if (targetPDBPath != "") {
        target = Structure(targetPDBPath);
        for (int i = 0; i < target.chainSize(); i++) targetChains.push_back(&target.getChain(i));
    }

    if (structureName != "") {
        Structure* seed = seedBin.getStructureNamed(structureName);
        for (Chain* C : targetChains) {
            Chain* Ccopy = new Chain(*C);
            seed->appendChain(Ccopy);
        }
        seed->writePDB(seed->getName()+".pdb");
        delete seed;
    }

    if (structureListPath != "") {
        vector<string> structureList = MstUtils::fileToArray(structureListPath);
        for (string structureName : structureList) {
            Structure* seed = seedBin.getStructureNamed(structureName);
            for (Chain* C : targetChains) {
                Chain* Ccopy = new Chain(*C);
                seed->appendChain(Ccopy);
            }
            seed->writePDB(seed->getName()+".pdb");
            delete seed;
        }
    }

    cout << "Done!" << endl;
    return 0;
}