#include "mstfasst.h"
#include "mstfuser.h"
#include "mstsystem.h"
#include "mstoptions.h"

#include "alignframes.h"
#include "chainextension.h"
#include "residuecontact.h"
#include "residueframe.h"
#include "scoreinterface.h"

int main(int argc, char *argv[]) {
    MstOptions op;
    op.setTitle("Extends either the N- or C-terminus of a provided binder fragment.");
    op.addOption("anchor", "The structure that will be extended",true);
    op.addOption("target", "The structure of the target protein",true);
    op.addOption("fasstDB", "Path to the FASST DB that will be searched to extend the chain",true);
    op.addOption("frameDB","Path to the frame DB that will be used to score chain extensions",true);
    op.addOption("potcontsJSON","Path to JSON file with 2D histograms describing potential contacts",true);
    op.addOption("terminus", "'N' or 'C' for N- or C-terminal extension",true);
    op.addOption("mode","The protocol to use to select the extending fragments: 'random', 'greedy', 'mc'. (default: greedy)");
    op.addOption("numOptions","The number of potential extensions that is sampled during each round when --mode = greedy (default: 10)",false);
    op.addOption("overlapSegmentLength","During each iteration, select this many residues from the terminus to define the extension query (default: 4)",false);
    op.addOption("extendSegmentLength", "During each iteration, extend the chain by this many residues (default: 4)",false);
    op.addOption("extendLength", "The number of residues to extend the structure by (default: 12)",false);
    op.setOptions(argc,argv);

    MstTimer timer;
    MstUtils::seedRandEngine(42);
    
    Structure anchor(op.getString("anchor"),"SKIPHETERO|ALLOW ILE CD1");
    Structure target(op.getString("target"),"SKIPHETERO|ALLOW ILE CD1");
    string structure_name = MstSys::splitPath(target.getName(),1);
    string fasstDBPath = op.getString("fasstDB");
    string frameDBPath = op.getString("frameDB");
    string potcontsJSONPath = op.getString("potcontsJSON");
    string terminus = op.getString("terminus");
    string mode = op.getString("mode","greedy");
    int numExtensions = op.getInt("numOptions",10);
    int overlap_length = op.getInt("overlapSegmentLength",4);
    int extend_stride = op.getInt("extendSegmentLength",4);
    int extend_length = op.getInt("extendLength",12);

    // extensionDirection extDir;
    // if ((terminus == "N")) {
    //     extDir = extensionDirection::NTERM;
    // } else if ((terminus == "C")) {
    //     extDir = extensionDirection::CTERM;
    // } else {
    //     MstUtils::error("Must specify either 'N' or 'C' for --terminus");
    // }
    // if ((mode != "random")&(mode != "greedy")&(mode != "mc")) MstUtils::error("Must specify either 'random', 'greedy', or 'mc' for --mode");

    // binderChainExtensionFASST chainExt(extDir,anchor,fasstDBPath);

    // binderScorerParams params;
    // params.frameDBPath = frameDBPath;
    // params.potentialContactsJSONPath = potcontsJSONPath;
    // params.posCut = 1.0;
    // params.oriCut = 10;
    // params.renormalizeProbabilities = false;

    // chainExt.setFixedStructuralContext(target,params);
    // chainExt.setOverlapSegmentLength(overlap_length);
    // chainExt.setExtensionSegmentLength(extend_stride);
    // if (mode == "random") {
    //     chainExt.extendChainRandom(extend_length);
    // } else if (mode == "greedy") {
    //     chainExt.extendChainGreedy(extend_length,numExtensions);
    // }
    return 0;
}