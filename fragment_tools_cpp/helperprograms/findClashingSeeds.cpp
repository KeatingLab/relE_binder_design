#include "mstoptions.h"
#include "mstrotlib.h"
#include "mstsystem.h"

#include "generateseeds.h"
#include "residuecontact.h"
#include "residueframe.h"
#include "utilities.h"
#include "utilitiesio.h"

int main(int argc, char *argv[]) {
    MstOptions op;
    op.setTitle("Loads a binary file and checks whether the given seed clashes with sidechain/backbone of another structure");
    op.addOption("seedBin","Path to seed binary file. In read mode, will load this file. In write mode will either A) load the file, if one already exists or B) create a new one",true);
    op.addOption("targetPDB","Path to a PDB structure to check for clashes",false);
    op.addOption("complexPDB","Path to a PDB structure to check for clashes",false);
    op.addOption("binderChainID","Chain ID of the binder in --complexPDB",false);
    op.setOptions(argc,argv);
        
    string seedBinPath = op.getString("seedBin");
    string targetPDBPath = op.getString("targetPDB","");
    string complexPDBPath = op.getString("complexPDB","");
    string binderChainID = op.getString("binderChainID","");

    if (!MstSys::fileExists(seedBinPath)) MstUtils::error("When in 'read' mode, must provide the path to an existing seed binary file path");

    seedBinaryFile seedBin(seedBinPath);
    seedBin.scanFilePositions();

    clashChecker checker;
    Structure target;
    if (targetPDBPath != "") {
        target = Structure(targetPDBPath);
        checker.setStructure(target);
    } else if ((complexPDBPath != "")&&(binderChainID != "")) {
        Structure complex = Structure(complexPDBPath);
        vector<Residue*> target_residues;
        cout << complex.residueSize() << " residues in complex" << endl;
        for (int i = 0; i < complex.chainSize(); i++) {
            Chain* C =  &complex.getChain(i);
            if (C->getID() == binderChainID) continue;
            vector<Residue*> chain_residues = C->getResidues();
            target_residues.insert(target_residues.end(),chain_residues.begin(),chain_residues.end());
        }
        target = Structure(target_residues);
        cout << target.residueSize() << " residues in target after deleting chain" << endl;
        checker.setStructure(target);
    } else {
        MstUtils::error("Must provide either --targetPDBPath or --complexPDBPath and --binderChainID");
    }

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