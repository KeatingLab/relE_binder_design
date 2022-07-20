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
    op.addOption("peptideRMSD","If provided, *instead* of storing to a binary file, will compute the backbone RMSD between each peptide structure and the reference (note: all peptides must have the same length)",false);
    op.setOptions(argc,argv);

    string pdbListPath = op.getString("pdbList");
    string binderChainID = op.getString("binderChainID");
    string outPrefix = op.getString("outPrefix");
    bool peptideRMSD = op.isGiven("peptideRMSD");

    string seedBinPath = outPrefix + "_binder.bin";

    fstream out;
    string fileName = outPrefix + "_info.csv";
    MstUtils::openFile(out,fileName,fstream::out);
    out << "outPrefix,binderName,binderNumRes,targetBBRMSD,binderBBRMSD" << endl;

    vector<string> pdbList = MstUtils::fileToArray(pdbListPath);
    if (pdbList.size() < 1) MstUtils::error("PDB list has less than one entry","main()");
    Structure firstS(pdbList[0]);
    cout << "First structure has " << firstS.chainSize() << " chains" << endl;
    for (int i = 0; i < firstS.chainSize(); i++) cout << "Chain ID: " << firstS.getChain(i) << endl;

    cout << "Isolating the target from the first structure: " << firstS.getName() << endl;
    vector<Atom*> refTargetAtoms = extractChainBBAtoms(firstS,binderChainID,false);
    vector<Atom*> refBinderAtoms = extractChainBBAtoms(firstS,binderChainID,true);
    Structure refTarget(refTargetAtoms);
    refTarget.writePDB(outPrefix+"_target.pdb");

    seedBinaryFile* seedBin = nullptr;
    if (!peptideRMSD) seedBin = new seedBinaryFile(seedBinPath,false);

    RMSDCalculator calc;
    for (string pdbPath : pdbList) {
        string sName = MstSys::splitPath(pdbPath,1);
        cout << "Reading " << sName << " ..." << endl;
        Structure S(pdbPath);
        vector<Atom*> allAtoms = S.getAtoms();
        vector<Atom*> targetAtoms = extractChainBBAtoms(S,binderChainID,false);
        vector<Atom*> binderAtoms = extractChainBBAtoms(S,binderChainID,true);

        // cout << "targetAtoms length: " << targetAtoms.size() << endl;
        // cout << "refTargetAtoms length: " << refTargetAtoms.size() << endl;
        // cout << "binderAtoms length: " << binderAtoms.size() << endl;

        calc.align(targetAtoms,refTargetAtoms,allAtoms);

        Structure transformedBinder(binderAtoms);
        // cout << "transformedBinder atom length: " << transformedBinder.atomSize() << endl;
        // cout << "transformedBinder residue length: " << transformedBinder.residueSize() << endl;
        transformedBinder.setName(sName);

        out << outPrefix << "," << sName << "," << transformedBinder.residueSize() << "," << calc.lastRMSD();

        if (seedBin != nullptr) {
            seedBin->appendStructure(&transformedBinder);
            out << "," << -1.0;
        } else {
            // Compute RMSD between peptide backbone atoms
            if (binderAtoms.size() != refBinderAtoms.size()) MstUtils::error("Mismatched number of backbone atoms between binder chain of "+pdbPath+" and the reference structure "+pdbList[0]);
            mstreal binderRMSD = calc.rmsd(binderAtoms,refBinderAtoms);
            out << "," << binderRMSD;
            S.writePDB(sName+"_aligned.pdb");
        }
        out << endl;
    }
    out.close();

    cout << "Done!" << endl;
    return 0;
}