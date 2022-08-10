#include "mstfuser.h"
#include "mstrotlib.h"
#include "msttypes.h"
#include "mstoptions.h"

#include "bridgeseeds.h"
#include "generateseeds.h"
#include "residueframe.h"
#include "utilities.h"

vector<int> getFragResIdx(int nTermRes, int resLength) {
    vector<int> result(resLength,0);
    iota(result.begin(),result.end(),nTermRes);
    return result;
}

int main (int argc, char *argv[]) {
    MstOptions op;
    op.setTitle("Given a list of seeds, reports all pairs of seeds that can be joined by a bridging-fragment.");
    op.addOption("seedBin","The path to a seed binary file",true);
    op.addOption("seedList","The path to the list of seeds that will be considered when searching for pairs of bridge-able seeds",false);
    op.addOption("bridgeDB","The path to a bridgeDB file. If provided, will load and search for matches to the terminus in the DB",true);
    op.addOption("structureDB","The path to a structureDB file. Required to verify that the matches are superimposable with the query",true);
    op.addOption("dCut","Cutoff applied when comparing CA distances (default: 0.5 Å");
    op.addOption("RMSDCut","Cutoff applied when calculating RMSD between superimposed termini (default 0.5 Å");
    op.addOption("selectedSeed","Either the path to PDB file or the name of a seed in the binary file. If provided, will ONLY search for bridges between this seed and the other seeds (either all seeds in the binary file, or just those in the list)",false);
    op.setOptions(argc,argv);

    string seedBinPath = op.getString("seedBin");
    string seedListPath = op.getString("seedList","");
    string bridgeDBPath = op.getString("bridgeDB");
    string structureDBPath = op.getString("structureDB");
    mstreal dCut = op.getReal("dCut",0.5);
    mstreal RMSDCut = op.getReal("RMSDCut",0.5);
    string selectedSeedName = op.getString("selectedSeed","");

    int overlapLength = 3;

    // Load the seeds from the binary file
    seedBinaryFile seeds(seedBinPath);
    seeds.scanFilePositions();

    // Load the bridge distances/structures from which the distances were calculated
    findSeedBridge bridgeFinder(op.getString("bridgeDB"));
    bridgeFinder.loadStructureDB(op.getString("structureDB"));
    bridgeFinder.reportAPVBoundaries();
    
    vector<Structure*> selectedSeeds;
    vector<string> seedList = MstUtils::fileToArray(seedListPath);
    for (string seedName : seedList) {
        Structure* selectedSeed = seeds.getStructureNamed(seedName);
        selectedSeeds.push_back(selectedSeed);
    }

    // For each pair of seeds, find whether they can be bridged in either direction
    fstream out;
    MstUtils::openFile(out,"bridgeSeeds.csv",fstream::out);
    out << "seedA,seedB,nRes,filename" << endl;

    string bridgeDirName = "seedBridgeStructures/";
    if (!MstSys::isDir(bridgeDirName)) MstSys::cmkdir(bridgeDirName);
    string fusedDirName = "fusedSeedBridgeStructures/";
    if (!MstSys::isDir(fusedDirName)) MstSys::cmkdir(fusedDirName);

    fusionParams params;
    fusionOutput fuserOut;

    MstTimer timer; timer.start();
    int distanceMatches = 0;
    int rmsdMatches = 0;
    string filename = "";
    for (int i = 0; i < selectedSeeds.size(); i++) {
        Structure* seedA = selectedSeeds[i];
        for (int j = 0; j < selectedSeeds.size(); j++) {
            Structure* seedB = selectedSeeds[j];

            bridgeFinder.setSearchQuery(seedA,seedB);
            timer.start();
            distanceMatches = bridgeFinder.searchSeedsByCADistance(dCut);
            timer.stop();
            cout << "Took " << timer.getDuration(MstTimer::msec) << " ms to find matches" << endl;
            // bridgeFinder.getBridgeLengthDist();
            if (distanceMatches == 0) continue;
            timer.start();
            rmsdMatches = bridgeFinder.verifyMatchesBySuperposition(RMSDCut);
            timer.stop();
            if (rmsdMatches == 0) continue;
            cout << "Took " << timer.getDuration(MstTimer::msec) << " ms to verify the matches" << endl;
            // bridgeFinder.getVerifiedBridgeLengthDist();
            vector<Structure> representativeBridges = bridgeFinder.getRepresentativeForEachLength();
            for (int k = 0; k < representativeBridges.size(); k++) {
                Structure& bridge = representativeBridges[k];
                int bridgeLength = bridge.residueSize() - (2*overlapLength);
                filename = "bridge_"+seedA->getName()+"_"+MstUtils::toString(bridgeLength)+"_"+seedB->getName();
                bridge.writePDB(bridgeDirName+filename+".pdb");
                out << seedA->getName() << "," << seedB->getName() << "," << bridgeLength << "," << filename << endl;

                // Fuse seeds + bridge and store
                fusionTopology topology(seedA->residueSize()+seedB->residueSize()+bridgeLength);
                topology.addFragment(*seedA,getFragResIdx(0,seedA->residueSize()));
                topology.addFragment(bridge,getFragResIdx(seedA->residueSize()-overlapLength,bridge.residueSize()));
                topology.addFragment(*seedB,getFragResIdx(seedA->residueSize()+bridgeLength,seedB->residueSize()));
                Structure fusedS = Fuser::fuse(topology,fuserOut,params);
                fusedS.writePDB(fusedDirName+filename+"_fused.pdb");

            }
        }
    }

    for (Structure* seed : selectedSeeds) delete seed;

    cout << "Done!" << endl;
}