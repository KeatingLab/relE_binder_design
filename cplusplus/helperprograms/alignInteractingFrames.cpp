#include "msttypes.h"
#include "mstoptions.h"

#include "alignframes.h"
#include "residuecontact.h"
#include "residueframe.h"
#include "utilities.h"

#include <chrono>

int main(int argc, char *argv[]) {
    MstOptions op;
    op.setTitle("Creates a mobile frame database from a protein structure database.");
    op.addOption("db","Path to protein structure database",false);
    op.addOption("pdb","Path to protein structure",false);
    op.addOption("resPair","If provided, will create a residue pair database instead of a mobile frame database",false);
    op.addOption("flankingRes","If provided, will also store pairs of residues within this many number of residues in the chain (default = -1)",false);
    op.addOption("flankingResSubsample","If provided, will subsample flanking residues at the provided rate (default = -1)",false);
    op.addOption("out","File name prefix",false);
    op.addOption("sample","Sample interacting residues at the specified rate",false);
    op.addOption("verbose","If provided, write the name of each pair of residues that was found to contact to stdout",false);
    op.setOptions(argc,argv);

    if (op.isGiven("db") == op.isGiven("pdb")) MstUtils::error("Must provide either '--db' or '--pdb', but not both");

    string proteinDBPath = op.getString("db");
    if (op.isGiven("pdb")) {
        Structure S(op.getString("pdb"));
        vdwContacts C(S.getResidues());
        augmentedStructureDB protDB;
        protDB.addTarget(S);
        protDB.setVDWContacts(0,C.getAllInteractingRes());
        proteinDBPath = "structure.db"; 
        protDB.writeDBFile(proteinDBPath);
    }

    string DBPath;
    if (!op.isGiven("resPair")) DBPath = op.getString("out") + ".frames.db";
    else DBPath = op.getString("out") + ".respairs.db";
    mstreal sampleRate = op.getReal("sample",0.0001);
    int flankingRes = op.getInt("flankingRes",-1);
    mstreal flankingResSubsampleRate = op.getReal("flankingResSubsampleRate",-1.0);

    alignInteractingFrames alignF(proteinDBPath,op.isGiven("verbose"),flankingRes,flankingResSubsampleRate);
    bool read = false;
    if (!op.isGiven("resPair")) {
        frameDB* frameBin = new frameDB(DBPath,read);

        set<string> aaTypes = SeqToolsExtension::getAANames();

        for (string aa : aaTypes) {
            cout << "amino acid: " << aa << endl;

            alignF.setAA(aa);
            alignF.findMobileFrames();
            cout << "protein structure DB has " << alignF.getNumInteracting() << " residue interactions from this residue type" << endl;

            alignF.writeAlignedInteractingResToPDB(aa+"_interactions.pdb",sampleRate);
            alignF.writeMobileFramesToBin(frameBin);
        }
        delete frameBin;
    } else {
        resPairDB* rpBin = new resPairDB(DBPath,read);

        set<string> aaTypes = SeqToolsExtension::getAANames();

        for (string aa : aaTypes) {
            cout << "amino acid: " << aa << endl;

            alignF.setAA(aa);
            alignF.findMobileFrames();
            cout << "protein structure DB has " << alignF.getNumInteracting() << " residue interactions from this residue type" << endl;

            alignF.writeAlignedInteractingResToPDB(aa+"_interactions.pdb",sampleRate);
            alignF.writeResiduePairsToBin(rpBin);
        }
        delete rpBin;
    }

    alignF.writeInteractionData(op.getString("out"));

    cout << "Done!" << endl;
    return 0;
}