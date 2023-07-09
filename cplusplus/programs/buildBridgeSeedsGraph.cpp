#include "mstfuser.h"
#include "mstrotlib.h"
#include "msttypes.h"
#include "mstoptions.h"

#include "bridgeseeds.h"
#include "generateseeds.h"
#include "residueframe.h"
// #include "utilities.h"
#include "utilitiesio.h"

vector<Residue*> getSeedResWithResection(Structure* seed, bool nterminus = true, int offset = 0) {
    vector<Residue*> all_res = seed->getResidues();
    vector<Residue*> ret_res;
    for (int i = 0; i < all_res.size(); i++) {
        if (nterminus && (i < offset)) {
            // skip n-terminal residue
        } else if (!nterminus && (i >= (seed->residueSize() - offset))) {
            // skip c-terminal residue
        } else {
            ret_res.push_back(all_res[i]);
        }
    }
    return ret_res;
}

int main (int argc, char *argv[]) {
//     MstOptions op;
//     op.setTitle("Given a list of seeds, reports all pairs of seeds that can be joined by a bridging-fragment.");
//     op.addOption("seedBin","The path to a seed binary file",false);
//     op.addOption("seedList","The path to the list of seeds that will be considered when searching for pairs of connectable seeds",false);
//     op.addOption("seedPDBsList","The path to a list of seed PDB files that will be considered when searching for pairs of connectable seeds");
//     op.addOption("bridgeDB","The path to a bridgeDB file. If provided, will load and search for matches to the terminus in the DB",true);
//     op.addOption("structureDB","The path to a structureDB file. Required to verify that the matches are superimposable with the query",true);
//     op.addOption("dCut","Cutoff applied when comparing CA distances (default: 0.5 Å");
//     op.addOption("RMSDCut","Cutoff applied when calculating RMSD between superimposed termini (default 0.5 Å");
//     op.addOption("resectLength","When searching seed termini, allow removal of terminal seed residues up to specified number. This allows for identification of bridges between non-terminal residues. (default = 0)",false);
//     op.addOption("minSeedLength","Seeds must retain at least this many residues (only applies when resectLength is > 0)");
//     // op.addOption("seedA","The path to a seed structure. If provided, will search for overlaps between this seed and other seeds OR seedB only, if that argument is provided",false);
//     // op.addOption("seedB","The path to a seed structure. If provided, will search for overlaps between this seed and other seeds OR seedA only, if that argument is provided",false);
//     op.addOption("seedANames","The path to a list of seed names in the binary file",false);
//     op.addOption("seedBNames","The path to a list of seed names in the binary file",false);
//     op.addOption("seedAPDBList","The path to a list of seed PDB files",false);
//     op.addOption("seedBPDBList","The path to a list of seed PDB files",false);
//     op.addOption("avoidClashesToStructure","Path to PDB file. If provided, will check for clash between briges and PDB and exclude bridges that clash",false);
//     op.addOption("ignoreSeedClash","",false);
//     op.addOption("maxSeedTerminiDistance","Seed termini will be searched for connections if they are within this distance (default: -1 Å)",false);
//     op.addOption("maxBridgeLength","Loop regions connecting seeds will contain no more than this many residues (default: 8 residues)",false);
//     op.addOption("workerID","The 0-indexed unique index for this job, e.g. ${SLURM_ARRAY_TASK_ID} (default = 0)",false);
//     op.addOption("nWorkers","The total number of jobs running in parallel (default = 1)",false);
//     op.addOption("verbose","");
//     // op.addOption("selectedSeed","Either the path to PDB file or the name of a seed in the binary file. If provided, will ONLY search for bridges between this seed and the other seeds (either all seeds in the binary file, or just those in the list)",false);
//     op.setOptions(argc,argv);

//     string seedBinPath = op.getString("seedBin","");
//     string seedListPath = op.getString("seedList","");
//     string seedPDBListPath = op.getString("seedPDBList","");
//     string bridgeDBPath = op.getString("bridgeDB");
//     string structureDBPath = op.getString("structureDB");
//     mstreal dCut = op.getReal("dCut",0.5);
//     mstreal RMSDCut = op.getReal("RMSDCut",0.5);
//     int resectLength = op.getInt("resectLength",0);
//     // string seedAPath = op.getString("seedA","");
//     // string seedBPath = op.getString("seedB","");
//     string seedANamesPath = op.getString("seedANames","");
//     string seedBNamesPath = op.getString("seedBNames","");
//     string seedAPDBListPath = op.getString("seedAPDBList","");
//     string seedBPDBListPath = op.getString("seedBPDBList","");
//     string clashStructure = op.getString("avoidClashesToStructure","");
//     bool ignoreSeedClash = op.isGiven("ignoreSeedClash");
//     int maxSeedTerminiDistance = op.getInt("maxSeedTerminiDistance",-1.0);
//     int maxBridgeLength = op.getInt("maxBridgeLength",8);
//     int workerID = op.getInt("workerID",0);
//     int nWorkers = op.getInt("nWorkers",1);
//     bool verbose = true; // op.isGiven("verbose");

//     int overlapLength = 3;
//     int minSeedLength = op.getInt("minSeedLength",overlapLength);

//     // if ((minimumSeedLength != 0)&(minimumSeedLength < overlapLength)) MstUtils::error("Minimum seed length :"+MstUtils::toString(minimumSeedLength)+" must be at least as great as overlap length: "+MstUtils::toString(overlapLength));
//     if ((seedListPath != "")&(seedBinPath=="")) MstUtils::error("If a seed list provided, a seed binary file must also be provided");

//     // // Load the seeds from the binary file
//     // seedBinaryFile* seeds = nullptr;
//     // if (seedBinPath != "") {
//     //     seeds = new seedBinaryFile(seedBinPath);
//     //     seeds->scanFilePositions();
//     // }

//     if (maxSeedTerminiDistance < 0) {
//         // Set the max distance according to simple rule:
//         // 3.3 Å / residue (distance between Ca in consecutive residues in extended conformation)
//         // Add two residues because the distance check is between central residues in the terminal regions
//         maxSeedTerminiDistance = 3.3 * (maxBridgeLength+2); 
//     }
    
//     // Load the bridge distances/structures from which the distances were calculated
//     findSeedBridge bridgeFinder(op.getString("bridgeDB"),MstUtils::toString(workerID),maxBridgeLength);
//     bridgeFinder.setMaxSeedDistance(maxSeedTerminiDistance);
//     bridgeFinder.loadStructureDB(op.getString("structureDB"));
//     bridgeFinder.reportAPVBoundaries();
    
//     seedPairDistributor seedPairs(workerID,nWorkers,true,seedBinPath);
//     if (((seedListPath!="")||(seedPDBListPath!=""))) {
//         // By default, perform all-to-all comparison
//         if (seedListPath != "") {
//         vector<string> seedList = MstUtils::fileToArray(seedListPath);
//         seedPairs.setSeedGroupSym(seedList);
//         }
//         if (seedPDBListPath != "") {
//             vector<string> seedPDBPathList = MstUtils::fileToArray(seedPDBListPath);
//             vector<shared_ptr<Structure>> seeds;
//             for (string path : seedPDBPathList) {
//                 shared_ptr<Structure> seed = make_shared<Structure>(path);
//                 seed->setName(structurePathToName(seed->getName()));
//                 seeds.push_back(seed);
//             }
//             seedPairs.setSeedGroupSym(seeds);
//         }
//     } else {
//         // In some cases, we will want to find connections between two distinct groups of seeds
//         vector<string> seed_names_A, seed_names_B;
//         vector<shared_ptr<Structure>> seeds_A, seeds_B;
//         if (seedAPDBListPath != "") {
//             vector<string> seedPDBPathList = MstUtils::fileToArray(seedAPDBListPath);
//             seeds_A = loadSharedPointerStructuresFromPaths(seedPDBPathList,"0");
//         }
//         if (seedBPDBListPath != "") {
//             vector<string> seedPDBPathList = MstUtils::fileToArray(seedBPDBListPath);
//             seeds_B = loadSharedPointerStructuresFromPaths(seedPDBPathList,"0");
//         }
//         if (seedANamesPath != "") {
//             seed_names_A = MstUtils::fileToArray(seedANamesPath);
//         } 
//         if (seedBNamesPath != "") {
//             seed_names_B = MstUtils::fileToArray(seedBNamesPath);
//         }
//         seedPairs.setSeedGroupAsym(seed_names_A,seed_names_B,seeds_A,seeds_B);
//     }

//     // string bridgeDirName = "seedBridgeStructures/";
//     // if (!MstSys::isDir(bridgeDirName)) MstSys::cmkdir(bridgeDirName);
//     // string fusedDirName = "fusedSeedBridgeStructures/";
//     // if (!MstSys::isDir(fusedDirName)) MstSys::cmkdir(fusedDirName);

//     MstTimer timer; timer.start();
//     int distanceMatches = 0;
//     int rmsdMatches = 0;
//     string filename = "";
//     string seed_names = "";
//     // int minSeedLengthA = minimumSeedLength;
//     // int minSeedLengthB = minimumSeedLength;

//     fuseSeedsAndBridge seedFuser(bridgeFinder.bridgeData.getTerminusLength(),MstUtils::toString(workerID));
//     if (clashStructure != "") {
//         Structure S(clashStructure);
//         seedFuser.setClashCheck(S);
//     }

//     clashChecker clashCheck;
//     cout << "Searching for connections between as many as " << seedPairs.batchSize() << " pairs of seeds" << endl;
//     int count = 0;
//     while (seedPairs.hasNext()) {
//         pair<Structure*,Structure*> seed_pair = seedPairs.next();
//         Structure* seedA = seed_pair.first;
//         Structure* seedB = seed_pair.second;
//         cout << count++ << " " << seedA->getName() << " and " << seedB->getName() << endl;
//         if (seedA == seedB) continue;
//         // Varying the offset allows us to "resect" the termini in case non-terminal residues are better for forming the bridge
//         for (int seedAOffset = 0; seedAOffset <= min(resectLength,max(0,seedA->residueSize()-minSeedLength)); seedAOffset++) {
//             if (verbose) cout << "seed A offset: " << seedAOffset << endl;
//             vector<Residue*> resectedSeedARes = getSeedResWithResection(seedA,false,seedAOffset);
//             if (!ignoreSeedClash) clashCheck.setResidues(resectedSeedARes);
//             for (int seedBOffset = 0; seedBOffset <= min(resectLength,max(0,seedB->residueSize()-minSeedLength)); seedBOffset++) {
//                 if (verbose) cout << "seed B offset: " << seedBOffset << endl;
                
//                 // Check for clash between seeds (with consideration of the offset)
//                 if (!ignoreSeedClash) {
//                     vector<Residue*> resectedSeedBRes = getSeedResWithResection(seedB,true,seedBOffset);
//                     if (clashCheck.checkForClashesToStructure(resectedSeedBRes)) {
//                         if (verbose) cout << "seeds clash" << endl;
//                         continue;
//                     }
//                 }

//                 // seed A -> seed B 
//                 bridgeFinder.setSearchQuery(seedA,seedB,seedAOffset,seedBOffset);
//                 timer.start();
//                 distanceMatches = bridgeFinder.searchSeedsByCADistance(dCut);
//                 timer.stop();
//                 if (verbose) std::cout << "Took " << timer.getDuration(MstTimer::msec) << " ms to find matches" << endl;
//                 if (distanceMatches == 0) continue;
//                 timer.start();
//                 rmsdMatches = bridgeFinder.verifyMatchesBySuperposition(RMSDCut);
//                 timer.stop();
//                 if (verbose) std::cout << "Took " << timer.getDuration(MstTimer::msec) << " ms to verify the matches" << endl;
//                 if (rmsdMatches == 0) continue;
//                 if (verbose) cout << "Found bridges between " << seedA->getName() << " with offset " << MstUtils::toString(seedAOffset) << " and " << seedB->getName() << " with offset " << MstUtils::toString(seedBOffset) << endl;

//                 // Save bridge structure info
//                 seed_names = seedA->getName()+"_"+MstUtils::toString(seedA->residueSize())+"_"+MstUtils::toString(seedAOffset)+"_"+MstUtils::toString(seedBOffset)+"_"+seedB->getName()+"_"+MstUtils::toString(seedB->residueSize());
//                 bridgeFinder.getVerifiedBridgeLengthDist(seed_names);
//                 // bridgeFinder.writeToMultiPDB(bridgeDirName+"/"+seed_names+".pdb",seed_names);

//                 // Fuse seeds with bridge
//                 seedFuser.setSeeds(seedA,seedB,seedAOffset,seedBOffset);
//                 seedFuser.setBridgeStructures(bridgeFinder.getVerifiedBridgeStructures());
//                 seedFuser.fuse();
//                 seedFuser.writeFusedStructuresToPDB();
//             }
//         }
//     }
//     cout << "Done!" << endl;
}