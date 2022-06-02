#include "mstfasst.h"
#include "mstfuser.h"
#include "mstsystem.h"
#include "mstoptions.h"

#include "alignframes.h"
#include "chainextension.h"
#include "fragmentdb.h"
#include "residuecontact.h"
#include "residueframe.h"
#include "scoreinterface.h"

int main(int argc, char *argv[]) {
    MstOptions op;
    op.setTitle("Extends either the N- or C-terminus of a provided binder fragment.");
    op.addOption("binder", "The structure of the binder chain",true);
    op.addOption("target", "The structure of the target protein",true);
    op.addOption("segmentOverlapGraph", "Path to the overlap graph with pre-computed termini overlap information",true);
    op.addOption("frameDB","Path to the frame DB that will be used to score chain extensions",true);
    op.addOption("potcontsJSON","Path to JSON file with 2D histograms describing potential contacts",true);
    op.addOption("terminus", "'N' or 'C' for N- or C-terminal extension",true);
    op.setOptions(argc,argv);

    MstTimer timer;
    MstUtils::seedRandEngine(42);
    
    Structure binder(op.getString("binder"),"SKIPHETERO|ALLOW ILE CD1");
    Structure target(op.getString("target"),"SKIPHETERO|ALLOW ILE CD1");
    string structure_name = MstSys::splitPath(target.getName(),1);
    string segmentOverlapGraphPath = op.getString("segmentOverlapGraph");
    string frameDBPath = op.getString("frameDB");
    string potcontsJSONPath = op.getString("potcontsJSON");
    string terminus = op.getString("terminus");

    extensionDirection extDir;
    if ((terminus == "N")) {
        extDir = extensionDirection::NTERM;
    } else if ((terminus == "C")) {
        extDir = extensionDirection::CTERM;
    } else {
        MstUtils::error("Must specify either 'N' or 'C' for --terminus");
    }

    binder.setName(MstSys::splitPath(binder.getName(),1));
    if (binder.chainSize() != 1) MstUtils::error("Binder structure must have a single chain","testExtendBinderChainSegmentDBCover::main");
    Chain* C = &binder.getChain(0);

    binderScorerParams params;
    params.frameDBPath = frameDBPath;
    params.potentialContactsJSONPath = potcontsJSONPath;
    params.posCut = 1.0;
    params.oriCut = 10;
    params.renormalizeProbabilities = false;

    binderChainExtension chainExt(target,segmentOverlapGraphPath,params);

    chainExt.coverChainWithExtensionSegments(C,extDir);

    return 0;
}