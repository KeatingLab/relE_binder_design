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
    op.addOption("anchor", "The structure that will be extended",true);
    op.addOption("target", "The structure of the target protein",true);
    op.addOption("segmentOverlapGraph", "Path to the overlap graph with pre-computed termini overlap information",true);
    op.addOption("frameDB","Path to the frame DB that will be used to score chain extensions",true);
    op.addOption("potcontsJSON","Path to JSON file with 2D histograms describing potential contacts",true);
    op.addOption("terminus", "'N' or 'C' for N- or C-terminal extension",true);
    op.addOption("extendLength","The number of residues to extend the structure by (default: 12)",false);
    op.addOption("mode","The protocol to use to select the extending fragments: 'random', 'greedy'. (default: greedy)");
    op.addOption("beamWidth","If in beam search mode, the number of extensions to select each round (default 5)",false);
    op.setOptions(argc,argv);

    MstTimer timer;
    MstUtils::seedRandEngine(42);
    
    Structure anchor(op.getString("anchor"),"SKIPHETERO|ALLOW ILE CD1");
    Structure target(op.getString("target"),"SKIPHETERO|ALLOW ILE CD1");
    string structure_name = MstSys::splitPath(target.getName(),1);
    string segmentOverlapGraphPath = op.getString("segmentOverlapGraph");
    string frameDBPath = op.getString("frameDB");
    string potcontsJSONPath = op.getString("potcontsJSON");
    string terminus = op.getString("terminus","");
    int16_t extendLength = op.getInt("extendLength",12);
    string mode = op.getString("mode","greedy");
    int beamWidth = op.getInt("beamWidth",5);

    // extensionDirection extDir;
    // if ((terminus == "N")) {
    //     extDir = extensionDirection::NTERM;
    // } else if ((terminus == "C")) {
    //     extDir = extensionDirection::CTERM;
    // } else {
    //     extDir = extensionDirection::EITHER;
    //     cout << "Will extend in both the N- and C-terminal directions" << endl;
    // }
    // if ((mode != "random")&(mode != "greedy")&(mode != "beam")) MstUtils::error("Must specify either 'random', 'greedy', or 'beam' for --mode");

    // binderScorerParams params;
    // params.frameDBPath = frameDBPath;
    // params.potentialContactsJSONPath = potcontsJSONPath;
    // params.posCut = 1.0;
    // params.oriCut = 10;
    // params.renormalizeProbabilities = false;

    // binderChainExtension chainExt(target,segmentOverlapGraphPath,params);
    // chainExt.setBinderAnchor(anchor);

    // if (mode == "random") {
    //     cout << "Extend mode: random" << endl;
    //     chainExt.extendChain(extendLength,extDir,extensionMode::RANDOM);
    // } else if (mode == "greedy") {
    //     cout << "Extend mode: greedy" << endl;
    //     chainExt.extendChain(extendLength,extDir,extensionMode::GREEDY);
    // } else if (mode == "beam") {
    //     cout << "Extend mode: beam" << endl;
    //     chainExt.extendChain(extendLength,extDir,extensionMode::BEAM,beamWidth);
    // }
    return 0;
}