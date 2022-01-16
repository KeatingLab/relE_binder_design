#include "msttypes.h"
#include "mstoptions.h"

#include "alignframes.h"
#include "scoreinterface.h"
#include "residuecontact.h"
#include "residueframe.h"

#include <chrono>

int main(int argc, char *argv[])
{
    MstOptions op;
    op.setTitle("Scores a protein-protein interface by finding the prevalence of pairwise interactions in the PDB");
    op.addOption("frameDB","Path to mobile frame database",true);
    op.addOption("structDB","Path to a database of protein structures",true);
    op.addOption("pdb","A PDB file defining the complex between the protein target and designed binder. The target must have a defined sequence",true);
    op.addOption("target","The protein chain(s) delimited by '_', e.g. 'A_B'",true);
    op.addOption("binder","The binder chain ID(s) delimited by '_', e.g. '0_1'",true);
    op.setOptions(argc,argv);

    string mobileFrameDB = op.getString("frameDB");
    string protDB = op.getString("structDB");
    string pdbPath = op.getString("pdb");
    string targetChainIDsString = op.getString("target");
    string binderChainIDsString = op.getString("binder");
    mstreal posCut = 0.5;
    mstreal oriCut = 15.0;

    // get chain ID lists
    vector<string> targetChainIDs = MstUtils::split(targetChainIDsString,"_");
    cout << "target chain IDs:";
    for (string s : targetChainIDs) cout << "\t" << s << endl;

    vector<string> binderChainIDs = MstUtils::split(binderChainIDsString,"_");
    cout << "binder chain IDs:";
    for (string s : binderChainIDs) cout << "\t" << s << endl;

    // initialize the scoring class
    interfaceScorer scorer(mobileFrameDB,posCut,oriCut);

    augmentedStructure complex(pdbPath);
    scorer.loadStructure(&complex,targetChainIDs,binderChainIDs);

    // score the interface
    scorer.scoreInterface();

    // write the score counts out to a file
    scorer.writeContactScoresToFile("score_test");

    // write the match structures
    alignInteractingFrames alignF(protDB);
    scorer.writeMatchStructures("match_structures.pdb",alignF);

    return 0;
}