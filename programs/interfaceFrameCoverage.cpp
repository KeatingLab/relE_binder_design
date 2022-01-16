#include "mstoptions.h"
#include "mstsystem.h"
#include "msttypes.h"

#include "alignframes.h"
#include "scoreinterface.h"
#include "residuecontact.h"
#include "residueframe.h"

#include <chrono>

int main(int argc, char *argv[])
{
    MstOptions op;
    op.setTitle("Defines protein-protein contacts and finds the prevalence of the pairwise interactions in the PDB");
    op.addOption("frameDB","Path to mobile frame database",true);
    op.addOption("structDB","Path to a database of protein structures",true);
    op.addOption("pdb","A PDB file defining the complex between the protein target and designed binder. The target must have a defined sequence",true);
    op.addOption("target","The protein chain ID(s) delimited by '_', e.g. 'A_B'",true);
    op.addOption("binder","The binder chain ID(s) delimited by '_', e.g. '0_1'",true);
    op.addOption("homThresh","The sequence homology threshold that will be applied when determining if a match is redundant. When the value is negative, no homology check is performed (default = 0.6)",false);
    op.setOptions(argc,argv);

    string mobileFrameDB = op.getString("frameDB");
    string protDB = op.getString("structDB","");
    string pdbPath = op.getString("pdb");
    string targetChainIDsString = op.getString("target");
    string binderChainIDsString = op.getString("binder");
    mstreal homThresh = op.getReal("homThresh",0.6);
    mstreal posCut = 0.5;
    mstreal oriCut = 15.0;

    // get chain ID lists
    vector<string> targetChainIDs = MstUtils::split(targetChainIDsString,"_");
    cout << "target chain IDs: ";
    for (string s : targetChainIDs) cout << " " << s;
    cout << endl;

    vector<string> binderChainIDs = MstUtils::split(binderChainIDsString,"_");
    cout << "binder chain IDs:";
    for (string s : binderChainIDs) cout << " " << s;
    cout << endl;

    // initialize the scoring class
    interfaceScorer scorer(mobileFrameDB,posCut,oriCut);
    alignInteractingFrames alignF(protDB);

    augmentedStructure complex(pdbPath);
    string complexName = MstSys::splitPath(complex.getName(),1);
    complex.setName(complexName);
    scorer.loadStructure(&complex,targetChainIDs,binderChainIDs);
    
    // score the interface
    cout << "Score the interface..." << endl;
    vdwContacts C(complex);
    C.setResidues(scorer.getTargetChains(),scorer.getBinderChains());
    vector<pair<Residue*,Residue*>> contacts = C.getInteractingRes();
    // for (auto pair : contacts) {
    //     cout << pair.first->getChainID() << pair.first->getNum() << " ";
    //     cout << pair.second->getChainID() << pair.second->getNum() << endl;
    // }
    if (homThresh > 0) {
        alignF.setHomThresh(homThresh);
        scorer.scoreInterface(contacts,&alignF);
    } else {
        scorer.scoreInterface(contacts);
    }

    // write the score counts out to a file
    cout << "Write score counts to file..." << endl;
    scorer.writeContactScoresToFile(complexName);
    scorer.writeContactMatchesToFile(complexName);
    scorer.writeContactPropertyToFile(complexName, interfaceScorer::property::COVERAGE);

    // write the match structures
    cout << "Write the match structures..." << endl;
    scorer.writeMatchStructures(complexName+"_match_structures.pdb",alignF);

    cout << "Done" << endl;
    return 0;
}