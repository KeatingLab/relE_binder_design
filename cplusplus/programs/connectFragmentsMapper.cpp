#include "mstfuser.h"
#include "mstrotlib.h"
#include "msttypes.h"
#include "mstoptions.h"

#include "bridgeseeds.h"
#include "generateseeds.h"
#include "residueframe.h"
#include "utilities.h"
#include "utilitiesio.h"

vector<Residue*> getSeedResWithResection(const shared_ptr<Structure>& seed, bool nterminus = true, int offset = 0) {
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
    MstOptions op;
    op.setTitle("Given a list of seeds, reports all pairs of seeds that can be joined by a bridging-fragment.");

    // symmetric mode (e.g. connecting a set of fragments to themselves)
    op.addOption("seedPDBsList","The path to a list of seed PDB files that will be considered when searching for pairs of connectable seeds");
    op.addOption("seedMultiPDB","The path to a seed multiPDB",false);
    op.addOption("seedList","The path to the list of seed names in multiPDB that will be considered when searching for pairs of connectable seeds",false);
    
    // asymmetric mode (e.g. connecting one set of fragments to another)
    op.addOption("seedAPDBList","The path to a list of seed PDB files",false);
    op.addOption("seedBPDBList","The path to a list of seed PDB files",false);
    op.addOption("seedAMultiPDB","The path to a seed multiPDB",false);
    op.addOption("seedBMultiPDB","The path to a seed multiPDB",false);
    op.addOption("seedANames","The path to a list of seed names in the multipdb file",false);
    op.addOption("seedBNames","The path to a list of seed names in the multipdb file",false);
    
    // other args/params
    op.addOption("bridgeDB","The path to a bridgeDB file. If provided, will load and search for matches to the terminus in the DB",true);
    op.addOption("structureDB","The path to a structureDB file. Required to verify that the matches are superimposable with the query",true);
    op.addOption("previousRunDir","The path to a directory where connectFragments was previously run. Provide this if you're connecting fragments that have already been fused and you want to preserve their topology information");
    op.addOption("dCut","Cutoff applied when comparing CA distances (default: 0.5 Å");
    op.addOption("RMSDCut","Cutoff applied when calculating RMSD between superimposed termini (default 0.5 Å");
    op.addOption("resectLength","When searching seed termini, allow removal of terminal seed residues up to specified number. This allows for identification of bridges between non-terminal residues. (default = 0)",false);
    op.addOption("minSeedLength","Seeds must retain at least this many residues (only applies when resectLength is > 0)");
    op.addOption("avoidClashesToStructure","Path to PDB file. If provided, will check for clash between briges and PDB and exclude bridges that clash",false);
    op.addOption("ignoreSeedClash","",false);
    op.addOption("maxSeedTerminiDistance","Seed termini will be searched for connections if they are within this distance (default: -1 Å)",false);
    op.addOption("maxBridgeLength","Loop regions connecting seeds will contain no more than this many residues (default: 8 residues)",false);
    op.addOption("workerID","The 0-indexed unique index for this job, e.g. ${SLURM_ARRAY_TASK_ID} (default = 0)",false);
    op.addOption("nWorkers","The total number of jobs running in parallel (default = 1)",false);
    op.addOption("verbose","");
    // op.addOption("selectedSeed","Either the path to PDB file or the name of a seed in the binary file. If provided, will ONLY search for bridges between this seed and the other seeds (either all seeds in the binary file, or just those in the list)",false);
    op.setOptions(argc,argv);

    MstTimer timer; timer.start();

    string seedPDBListPath = op.getString("seedPDBList","");
    string seedMultiPDBPath = op.getString("seedMultiPDB","");
    string seedListPath = op.getString("seedList","");

    string seedAPDBListPath = op.getString("seedAPDBList","");
    string seedBPDBListPath = op.getString("seedBPDBList","");
    string seedAMultiPDBPath = op.getString("seedAMultiPDB","");
    string seedBMultiPDBPath = op.getString("seedBMultiPDB","");
    string seedANamesPath = op.getString("seedANames","");
    string seedBNamesPath = op.getString("seedBNames","");

    string bridgeDBPath = op.getString("bridgeDB");
    string structureDBPath = op.getString("structureDB");
    string previousRunDirPath = op.getString("previousRunDir","");
    mstreal dCut = op.getReal("dCut",0.5);
    mstreal RMSDCut = op.getReal("RMSDCut",0.5);
    int resectLength = op.getInt("resectLength",0);
    string clashStructure = op.getString("avoidClashesToStructure","");
    bool ignoreSeedClash = op.isGiven("ignoreSeedClash");
    int maxSeedTerminiDistance = op.getInt("maxSeedTerminiDistance",-1.0);
    int maxBridgeLength = op.getInt("maxBridgeLength",8);
    int workerID = op.getInt("workerID",0);
    int nWorkers = op.getInt("nWorkers",1);
    bool verbose = true; // op.isGiven("verbose");

    int overlapLength = 3;
    int minSeedLength = op.getInt("minSeedLength",overlapLength);

    //// check that inputs are consistent with a single mode
    bool sym_args = ((seedPDBListPath!="")||(seedMultiPDBPath!="")||(seedListPath!=""));
    bool asym_args = ((seedAPDBListPath!="")||(seedBPDBListPath!="")||(seedAMultiPDBPath!="")||(seedBMultiPDBPath!="")||(seedANamesPath!="")||(seedBNamesPath!=""));
    if (sym_args==asym_args) MstUtils::error("Must provide either symmetric or asymmetric arguments, but not both:"+MstUtils::toString(sym_args)+","+MstUtils::toString(asym_args));
    if (sym_args&&(((seedListPath!="")||(seedMultiPDBPath!=""))==(seedPDBListPath!=""))) MstUtils::error("In sym mode, provide either multipdb + seed names OR list of paths to structures");
    if (asym_args&&(((seedANamesPath!="")||(seedAMultiPDBPath!=""))==(seedAPDBListPath!=""))) MstUtils::error("In asym mode, provide either multipdb + seed names OR list of paths to structures for group A");
    if (asym_args&&(((seedBNamesPath!="")||(seedBMultiPDBPath!=""))==(seedBPDBListPath!=""))) MstUtils::error("In asym mode, provide either multipdb + seed names OR list of paths to structures for group B");
    if (maxSeedTerminiDistance < 0) {
        // Set the max distance according to simple rule:
        // 3.3 Å / residue (distance between Ca in consecutive residues in extended conformation)
        // Add two residues because the distance check is between central residues in the terminal regions
        maxSeedTerminiDistance = 3.3 * (maxBridgeLength+2); 
    }
    
    //// Load DBs, seeds, and topology data
    // Load the bridge distances/structures from which the distances were calculated
    findSeedBridge bridgeFinder(op.getString("bridgeDB"),MstUtils::toString(workerID),maxBridgeLength);
    bridgeFinder.setMaxSeedDistance(maxSeedTerminiDistance);
    bridgeFinder.loadStructureDB(op.getString("structureDB"));
    bridgeFinder.reportAPVBoundaries();
    
    seedPairDistributor seedPairs(workerID,nWorkers,true);
    if (sym_args) {
        multiPDBFile seedMulti(seedMultiPDBPath);
        vector<shared_ptr<Structure>> seeds;
        if (seedListPath != "") {
            // Load from multiPDB
            vector<string> seedList = MstUtils::fileToArray(seedListPath);
            seeds = seedMulti.getStructuresSP(seedList);
        } else {
            // Load from paths
            vector<string> seedPDBPathList = MstUtils::fileToArray(seedPDBListPath);
            seeds = loadSharedPointerStructuresFromPaths(seedPDBPathList,"0");
        }
        seedPairs.setSeedGroupSym(seeds);
    } else {
        // In some cases, we will want to find connections between two distinct groups of seeds
        multiPDBFile seedAMulti(seedAMultiPDBPath), seedBMulti(seedBMultiPDBPath);
        vector<shared_ptr<Structure>> seeds_A, seeds_B;
        if (seedAPDBListPath != "") {
            vector<string> seedPDBPathList = MstUtils::fileToArray(seedAPDBListPath);
            seeds_A = loadSharedPointerStructuresFromPaths(seedPDBPathList,"0");
        } else {
            if (seedANamesPath == "") {
                seeds_A = seedAMulti.loadAllStructuresSP();
            } else {
                vector<string> seed_A_names = MstUtils::fileToArray(seedANamesPath);
                seeds_A = seedAMulti.getStructuresSP(seed_A_names);
            }
        }

        if (seedBPDBListPath != "") {
            vector<string> seedPDBPathList = MstUtils::fileToArray(seedBPDBListPath);
            seeds_B = loadSharedPointerStructuresFromPaths(seedPDBPathList,"0");
        } else {
            if (seedBNamesPath == "") {
                seeds_B = seedBMulti.loadAllStructuresSP();
            } else {
                vector<string> seed_B_names = MstUtils::fileToArray(seedBNamesPath);
                seeds_B = seedBMulti.getStructuresSP(seed_B_names);
            }
        }
        seedPairs.setSeedGroupAsym(seeds_A,seeds_B);
    }

    //// Create database files
    string fragmentDBname = "fragmentDB_"+MstUtils::toString(workerID)+".pdb";
    string fusedDBname = "fusedDB_"+MstUtils::toString(workerID)+".pdb";
    multiPDBFile fragmentsDB(fragmentDBname,false);
    multiPDBFile fusedDB(fusedDBname,false);
    string topoDBname = "topoDB_"+MstUtils::toString(workerID)+".json";
    topologyDB topoDB;

    // If some seeds are derived from previous fusion topologies, we load them now 
    map<string,shared_ptr<fragmentTopology>> previousTopologies;
    if (previousRunDirPath != "") {
        // load old topoDB
        string oldTopoDBpath = previousRunDirPath + "/topoDB.json";
        topologyDB oldTopoDB;
        oldTopoDB.readDBFromJSON(oldTopoDBpath);

        // load old fragDB
        string oldFragDBpath = previousRunDirPath + "/fragmentDB.pdb";
        multiPDBFile oldFragDB(oldFragDBpath);
        oldFragDB.countStructures();

        // check if any of the 'seeds' (which correspond to fusion topologies) are in the topoDB
        vector<string> all_seed_names = seedPairs.getAllSeedNames();
        for (string seed_name : all_seed_names) {
            if (oldTopoDB.isTopoInDB(seed_name)) {
                // if in DB, copy the fragments, create fragmentTopology, add to map
                shared_ptr<fragmentTopology> fT = oldTopoDB.getTopologyFromDB(seed_name,oldFragDB);
                previousTopologies[seed_name] = fT;
            }
        }

    }

    // Create info file for this job
    if (workerID == 0) {
        // Write basic information about the job to an extra file
        fstream job_info;
        MstUtils::openFile(job_info,"job_info.txt",fstream::out);
        job_info << "multipdb " << seedMultiPDBPath << endl;
        job_info << "multipdbA " << seedAMultiPDBPath << endl;
        job_info << "multipdbB " << seedBMultiPDBPath << endl;
        job_info << "structuredb " << structureDBPath << endl;
        job_info << "bridgedb" << bridgeDBPath << endl;
        job_info << "nworkers " << nWorkers << endl;
        job_info.close();
    }

    timer.stop();
    cout << "Startup took " << timer.getDuration() << " s" << endl;
     
    //// Connect fragments
    int distanceMatches = 0;
    int rmsdMatches = 0;
    string filename = "";
    string seed_names = "";

    clashChecker seedClashCheck, targetClashCheck;
    Structure targetS;
    if (clashStructure != "") {
        targetS = Structure(clashStructure);
        targetClashCheck.setStructure(targetS);
    }

    cout << "Searching for connections between as many as " << seedPairs.batchSize() << " pairs of seeds" << endl;
    int count = 0;
    while (seedPairs.hasNext()) {
        auto seed_pair = seedPairs.next();
        shared_ptr<Structure> seedA = seed_pair.first;
        shared_ptr<Structure> seedB = seed_pair.second;
        cout << "Progress: " << count++ << " " << seedA->getName() << " and " << seedB->getName() << endl;
        if (seedA == seedB) continue;
        // Varying the offset allows us to "resect" the termini in case non-terminal residues are better for forming the bridge
        for (int seedAOffset = 0; seedAOffset <= min(resectLength,max(0,seedA->residueSize()-minSeedLength)); seedAOffset++) {
            if (verbose) cout << "seed A offset: " << seedAOffset << endl;
            vector<Residue*> resectedSeedARes = getSeedResWithResection(seedA,false,seedAOffset);
            if (!ignoreSeedClash) seedClashCheck.setResidues(resectedSeedARes);
            for (int seedBOffset = 0; seedBOffset <= min(resectLength,max(0,seedB->residueSize()-minSeedLength)); seedBOffset++) {
                if (verbose) cout << "seed B offset: " << seedBOffset << endl;
                
                // Check for clash between seeds (with consideration of the offset)
                if (!ignoreSeedClash) {
                    vector<Residue*> resectedSeedBRes = getSeedResWithResection(seedB,true,seedBOffset);
                    if (seedClashCheck.checkForClashesToStructure(resectedSeedBRes)) {
                        if (verbose) cout << "seeds clash" << endl;
                        continue;
                    }
                }

                // seed A -> seed B 
                bridgeFinder.setSearchQuery(seedA,seedB,seedAOffset,seedBOffset);
                timer.start();
                distanceMatches = bridgeFinder.searchSeedsByCADistance(dCut);
                timer.stop();
                if (verbose) std::cout << "Took " << timer.getDuration(MstTimer::msec) << " ms to find matches" << endl;
                if (distanceMatches == 0) continue;
                timer.start();
                rmsdMatches = bridgeFinder.verifyMatchesBySuperposition(RMSDCut);
                timer.stop();
                if (verbose) std::cout << "Took " << timer.getDuration(MstTimer::msec) << " ms to verify the matches" << endl;
                if (rmsdMatches == 0) continue;
                if (verbose) cout << "Found bridges between " << seedA->getName() << " with offset " << MstUtils::toString(seedAOffset) << " and " << seedB->getName() << " with offset " << MstUtils::toString(seedBOffset) << endl;

                // Save bridge structure info
                seed_names = seedA->getName()+"_"+MstUtils::toString(seedA->residueSize())+"_"+MstUtils::toString(seedAOffset)+"_"+MstUtils::toString(seedBOffset)+"_"+seedB->getName()+"_"+MstUtils::toString(seedB->residueSize());
                bridgeFinder.getVerifiedBridgeLengthDist(seed_names);
                // bridgeFinder.writeToMultiPDB(bridgeDirName+"/"+seed_names+".pdb",seed_names);

                // Fuse seeds with bridge
                vector<vector<shared_ptr<Structure>>> clustered_seed_bridges_by_length = bridgeFinder.getVerifiedClusteredBridgeStructures();
                for (int bridge_length = 0; bridge_length < clustered_seed_bridges_by_length.size(); bridge_length++) {
                    for (int bridgeN = 0; bridgeN < clustered_seed_bridges_by_length[bridge_length].size(); bridgeN++) {
                        shared_ptr<Structure> connector = clustered_seed_bridges_by_length[bridge_length][bridgeN];
                        // Check for clashes to target
                        if ((targetClashCheck.isStructureSet())&&(targetClashCheck.checkForClashesToStructure(connector->getResidues()))) {
                            cout << "Bridge " << bridgeN << " clashes and will be discarded" << endl;
                            continue;
                        }
                        // Fuse the seeds and save to file
                        fragmentTopology topo;

                        // Add seed A (or the entire topology from which it was derived)
                        if (previousTopologies.count(seedA->getName()) == 0) {
                            topo.addInitialFragment(seedA,3,0,seedAOffset);
                        } else {
                            topo.addMultipleTopologyElements(previousTopologies[seedA->getName()]->getFragments());
                        }
                        // Add connector 
                        topo.addFragmentToEnd(connector,3,true,0,0,false);
                        // Add seed B (or the entire topology from which it was derived)
                        if (previousTopologies.count(seedB->getName()) == 0) {
                            topo.addFragmentToEnd(seedB,3,true,seedBOffset,0);
                        } else {
                            topo.addMultipleTopologyElements(previousTopologies[seedB->getName()]->getFragments());
                        }
                        topo.reportTopology();
                        shared_ptr<Structure> fused = topo.fuseTopology();

                        // Write fused structure and fragments to multiPDBfiles
                        fusedDB.addStructure(*fused,true); // each fused structure should be unique
                        topo.addFragmentsInTopoToDB(fragmentsDB);

                        // Add topology metadata to DB that will later be written
                        topoDB.addTopologyToDB(topo);
                    }
                }
            }
        }
    }
    topoDB.writeDBtoJSON(topoDBname);
    cout << "Done!" << endl;
}