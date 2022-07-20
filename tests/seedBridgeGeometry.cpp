#include "mstrotlib.h"
#include "msttypes.h"
#include "mstoptions.h"

#include "bridgeseeds.h"
#include "residueframe.h"
#include "utilities.h"

int main (int argc, char *argv[]) {
    MstOptions op;
    op.setTitle("Reports geometric parameters that reflect whether two segments can be easily joined by building in intervening residues");
    op.addOption("fragmentA","The PDB file describing the first fragment",true);
    op.addOption("fragmentB","The PDB file describing the second fragment",true);
    op.addOption("bridgeDB","The path to a bridgeDB file. If provided, will load and search for matches to the terminus in the DB",false);
    op.addOption("structureDB","The path to a structureDB file. Required to verify that the matches are superimposable with the query",false);
    op.addOption("dCut","Cutoff applied when comparing CA distances");
    op.addOption("RMSDCut","Cutoff applied when calculating RMSD between superimposed termini");
    op.setOptions(argc,argv);

    string fragmentApdbPath = op.getString("fragmentA");
    string fragmentBpdbPath = op.getString("fragmentB");

    Structure fragmentA(fragmentApdbPath);
    Structure fragmentB(fragmentBpdbPath);

    // get the terminal residue of each fragment
    Residue* fAntermR = fragmentA.getResidues().front();
    Residue* fActermR = fragmentA.getResidues().back();
    Residue* fBntermR = fragmentB.getResidues().front();
    Residue* fBctermR = fragmentB.getResidues().back();

    // fA C-terminus -> fB N-terminus
    cout << "fA C-terminus -> fB N-terminus" << endl;

    mstreal distance = seedResidueBridgeGeometry::residueDistance(fActermR,fBntermR);
    mstreal ctermAngleDeviation = seedResidueBridgeGeometry::cTerminalResidueDeviation(fActermR,fBntermR);
    mstreal ntermAngleDeviation = seedResidueBridgeGeometry::nTerminalResidueDeviation(fActermR,fBntermR);

    cout << "distance: " <<  distance << endl;
    cout << "ctermCosAngleDeviation: " << ctermAngleDeviation << endl;
    cout << "ntermCosAngleDeviation: " << ntermAngleDeviation << endl;
    cout << "ctermAngleDeviation: " << acos(ctermAngleDeviation)*(180/M_PI) << endl;
    cout << "ntermAngleDeviation: " << acos(ntermAngleDeviation)*(180/M_PI) << endl;


    // fB C-terminus -> fA C-terminus
    cout << "fB C-terminus -> fA C-terminus" << endl;

    distance = seedResidueBridgeGeometry::residueDistance(fBctermR,fAntermR);
    ctermAngleDeviation = seedResidueBridgeGeometry::cTerminalResidueDeviation(fBctermR,fAntermR);
    ntermAngleDeviation = seedResidueBridgeGeometry::nTerminalResidueDeviation(fBctermR,fAntermR);

    cout << "distance: " <<  distance << endl;
    cout << "ctermCosAngleDeviation: " << ctermAngleDeviation << endl;
    cout << "ntermCosAngleDeviation: " << ntermAngleDeviation << endl;
    cout << "ctermAngleDeviation: " << acos(ctermAngleDeviation)*(180/M_PI) << endl;
    cout << "ntermAngleDeviation: " << acos(ntermAngleDeviation)*(180/M_PI) << endl;

    if (op.isGiven("bridgeDB")) {
        findSeedBridge bridgeFinder(op.getString("bridgeDB"));
        bridgeFinder.loadStructureDB(op.getString("structureDB"));
        bridgeFinder.reportAPVBoundaries();
        bridgeFinder.setSearchQuery(&fragmentA,&fragmentB);
        MstTimer timer; timer.start();
        bridgeFinder.searchSeedsByCADistance(op.getReal("dCut",0.25));
        timer.stop();
        cout << "Took " << timer.getDuration(MstTimer::msec) << " ms to find matches" << endl;
        bridgeFinder.getBridgeLengthDist();
        timer.start();
        bridgeFinder.verifyMatchesBySuperposition(op.getReal("RMSDCut",0.25));
        timer.stop();
        cout << "Took " << timer.getDuration(MstTimer::msec) << " ms to verify the matches" << endl;
        bridgeFinder.getVerifiedBridgeLengthDist();
        vector<Structure> verifiedSeedBridges = bridgeFinder.getVerifiedBridgeStructures();
        cout << verifiedSeedBridges.size() << " v matches" << endl;
        if (!verifiedSeedBridges.empty()) verifiedSeedBridges.front().writePDB("seedBridge.pdb");
    }
    cout << "Done!" << endl;
}