#include "msttypes.h"
#include "mstoptions.h"

#include "alignframes.h"
#include "residuecontact.h"
#include "residueframe.h"

#include <chrono>

int main(int argc, char *argv[])
{
    MstOptions op;
    op.setTitle("Creates a mobile frame database from a protein structure database.");
    op.addOption("db","Path to protein structure database",false);
    op.addOption("pdb","Path to protein structure",false);
    op.addOption("out","File name prefix",false);
    op.addOption("sample","Sample interacting residues at the specified rate",false);
    op.setOptions(argc,argv);

    if (op.isGiven("db") == op.isGiven("pdb")) MstUtils::error("Must provide either '--db' or '--pdb', but not both");

    string proteinDBPath = op.getString("db");
    if (op.isGiven("pdb")) {
        Structure S(op.getString("pdb"));
        vdwContacts C(S);
        proteinFrameDB protDB;
        protDB.addTarget(S);
        protDB.setVDWContacts(0,C.getAllInteractingRes());
        proteinDBPath = "structure.db"; 
        protDB.writeDBFile(proteinDBPath);
    }

    string mobileFrameDBPath = op.getString("out") + ".frames.db";
    mstreal sampleRate = op.getReal("sample",0.0001);

    alignInteractingFrames alignF(proteinDBPath);
    bool read = false;
    frameDB* frameBin = new frameDB(mobileFrameDBPath,read);

    vector<string> aaTypes = {"ALA","ARG","ASN","ASP",
                             "CYS","GLN","GLU","GLY",
                             "HIS","ILE","LEU","LYS",
                             "MET","PHE","PRO","SER",
                             "THR","TRP","TYR","VAL"};

    for (string aa : aaTypes) {
        cout << "amino acid: " << aa << endl;

        alignF.setAA(aa);
        alignF.findMobileFrames();
        cout << "protein structure DB has " << alignF.getNumInteracting() << " residue interactions from this residue type" << endl;

        alignF.writeAlignedInteractingResToPDB(aa+"_interactions.pdb",sampleRate);
        alignF.writeMobileFramesToBin(frameBin);
    }
    alignF.writeInteractionData(op.getString("out"));

    delete frameBin;

    return 0;
}