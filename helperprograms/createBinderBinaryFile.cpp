#include "mstoptions.h"
#include "mstrotlib.h"
#include "mstsystem.h"

#include "generateseeds.h"
#include "residuecontact.h"
#include "residueframe.h"
#include "utilities.h"

vector<Atom*> extractChainBBAtoms(Structure& S, string chainID, bool binder) {
    vector<Chain*> selectedChains;
        Chain* binderC = S.getChainByID(chainID);
        if (binderC == NULL) MstUtils::error("Binder chain "+chainID+" not found in structure "+S.getName(),"extractChainBBAtoms()");
    if (binder) {
        selectedChains.push_back(binderC);
    } else {
        for (int i = 0; i < S.chainSize(); i++) {
            Chain* selectedChain = &S.getChain(i);
            if (selectedChain->getID() != chainID) selectedChains.push_back(selectedChain);
        }
    }
    return MiscTools::getBackboneAtoms(selectedChains);
}

int main(int argc, char *argv[]) {
    MstOptions op;
    op.setTitle("Adds designed binders to a binary file that can be scored downstream");
    op.addOption("pdbList","List of PDB structures",true);
    op.addOption("binderChainID","The chain ID of the binder",true);
    op.addOption("outPrefix","A prefix for naming all files",true);
    op.setOptions(argc,argv);

    string pdbListPath = op.getString("pdbList");
    string binderChainID = op.getString("binderChainID");
    string outPrefix = op.getString("outPrefix");

    string seedBinPath = outPrefix + "_binder.bin";
    seedBinaryFile seedBin(seedBinPath,false);

    fstream out;
    string fileName = outPrefix + "_info.csv";
    MstUtils::openFile(out,fileName,fstream::out);
    out << "outPrefix,binderName,binderNumRes,targetBBRMSD" << endl;

    vector<string> pdbList = MstUtils::fileToArray(pdbListPath);
    if (pdbList.size() < 1) MstUtils::error("PDB list has less than one entry","main()");
    Structure firstS(pdbList[0]);
    cout << "First structure has " << firstS.chainSize() << " chains" << endl;
    for (int i = 0; i < firstS.chainSize(); i++) cout << "Chain ID: " << firstS.getChain(i) << endl;

    cout << "Isolating the target from the first structure: " << firstS.getName() << endl;
    vector<Atom*> refTargetAtoms = extractChainBBAtoms(firstS,binderChainID,false);
    Structure refTarget(refTargetAtoms);
    refTarget.writePDB(outPrefix+"_target.pdb");

    RMSDCalculator calc;
    for (string pdbPath : pdbList) {
        string sName = MstSys::splitPath(pdbPath,1);
        cout << "Reading " << sName << " ..." << endl;
        Structure S(pdbPath);
        vector<Atom*> targetAtoms = extractChainBBAtoms(S,binderChainID,false);
        vector<Atom*> binderAtoms = extractChainBBAtoms(S,binderChainID,true);

        // cout << "targetAtoms length: " << targetAtoms.size() << endl;
        // cout << "refTargetAtoms length: " << refTargetAtoms.size() << endl;
        // cout << "binderAtoms length: " << binderAtoms.size() << endl;

        calc.align(targetAtoms,refTargetAtoms,binderAtoms);

        Structure transformedBinder(binderAtoms);
        // cout << "transformedBinder atom length: " << transformedBinder.atomSize() << endl;
        // cout << "transformedBinder residue length: " << transformedBinder.residueSize() << endl;
        transformedBinder.setName(sName);
        seedBin.appendStructure(&transformedBinder);

        out << outPrefix << "," << sName << "," << transformedBinder.residueSize() << "," << calc.lastRMSD() << endl;
    }
    out.close();

    cout << "Done!" << endl;
    return 0;
}