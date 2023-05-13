#include "mstfuser.h"
#include "mstrotlib.h"
#include "msttypes.h"
#include "mstoptions.h"

#include "bridgeseeds.h"
#include "generateseeds.h"
#include "residueframe.h"
#include "utilities.h"

class seedPairDistributor {
    // Note: object takes ownership of all seeds that are passed to it
    public:
        seedPairDistributor(int _workerID = 0, int _nWorkers = 1, string seedbinpath = "") {
            workerID = _workerID;
            nWorkers = _nWorkers;
            if (seedbinpath != "") {
                seedbin = new seedBinaryFile(seedbinpath);
                seedbin->scanFilePositions();
                cout << seedbin->structureCount() << " structures in binary file" << endl;
            }
        }

        ~seedPairDistributor() {
            for (Structure* seed: all_loaded_seeds) delete seed;
            if (seedbin != nullptr) delete seedbin;
        };

        // Directly provide the seed structures with the following three methods
        void setSeeds(vector<Structure*> seeds) {
            seed_group_a = seeds;
            seed_group_b = seeds;
            for (Structure* seed: seeds) all_loaded_seeds.insert(seed);
        };

        void setSeedsGroupA(vector<Structure*> seeds) {
            seed_group_a = seeds;
            for (Structure* seed: seeds) all_loaded_seeds.insert(seed);
        }

        void setSeedsGroupB(vector<Structure*> seeds) {
            seed_group_b = seeds;
            for (Structure* seed: seeds) all_loaded_seeds.insert(seed);
        }

        // Load the seeds from a binary file with the following methods
        void setSeeds(vector<string> seed_names) {
            vector<Structure*> seeds = seedbin->getStructuresNamed(seed_names);
            setSeeds(seeds);
        };

        void setSeedsGroupA(vector<string> seed_names) {
            vector<Structure*> seeds = seedbin->getStructuresNamed(seed_names);
            setSeedsGroupA(seeds);
        }

        void setSeedsGroupB(vector<string> seed_names) {
            vector<Structure*> seeds = seedbin->getStructuresNamed(seed_names);
            setSeedsGroupB(seeds);
        }

    bool hasNext() {
        // Batch is done when we've iterated over all 
        return ((group_a_current_index<group_a_max_index)|(group_b_current_index<group_b_max_index));
    }

    pair<Structure*,Structure*> next() {
        if ((group_a_current_index == -1)|(group_b_current_index == -1)) {
            // we round up when determining the batch size, so the final batch is smaller

            // There are A x B jobs total (imagine a grid). We use the worker ID and batch size to get our index in the total number of 
            // jobs (which is some value less than A * B - 1). We convert this into a position in the grid, e.g., we get the row/column
            // of A and B respectively
            int batch_size = batchSize();
            group_a_current_index = int(workerID * batch_size / seed_group_b.size());
            group_b_current_index = int(workerID * batch_size % seed_group_b.size());
            // group_a_min_index = group_a_current_index;
            // group_b_min_index = group_b_current_index;
            group_a_max_index = min(int(((workerID + 1) * batch_size) / seed_group_b.size()),int(seed_group_a.size()));
            group_b_max_index = min(int(((workerID + 1) * batch_size) % seed_group_b.size()),int(seed_group_b.size()));
        }
        pair<Structure*,Structure*> result(seed_group_a[group_a_current_index],seed_group_b[group_b_current_index]);
        // Now increment the position in the grid
        group_b_current_index++;
        if (group_b_current_index == seed_group_b.size()) {
            group_b_current_index = 0;
            group_a_current_index++;
        }
        return result;
    }

    int batchSize() {
        return ceil(mstreal(seed_group_a.size() * seed_group_b.size()) / mstreal(nWorkers));
    }

    private:
        int workerID = 0; // workers are 0-indexed
        int nWorkers = 1;

        int group_a_current_index = -1;
        int group_b_current_index = -1;

        // range is [min,max), size of range depends on number of seeds in each group and number of workers
        // idx is 0-indexed of course (to match the vector)
        int group_a_min_index = 0;
        int group_a_max_index = 0;
        int group_b_min_index = 0;
        int group_b_max_index = 0;

        seedBinaryFile* seedbin = nullptr;
        vector<Structure*> seed_group_a;
        vector<Structure*> seed_group_b;
        set<Structure*> all_loaded_seeds;
        
};

string seedPathToName(string seed) {
    return MstSys::splitPath(seed,1);
}

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
    MstOptions op;
    op.setTitle("Given a list of seeds, reports all pairs of seeds that can be joined by a bridging-fragment.");
    op.addOption("seedBin","The path to a seed binary file",false);
    op.addOption("seedList","The path to the list of seeds that will be considered when searching for pairs of connectable seeds",false);
    op.addOption("seedPDBsList","The path to a list of seed PDB files that will be considered when searching for pairs of connectable seeds");
    op.addOption("bridgeDB","The path to a bridgeDB file. If provided, will load and search for matches to the terminus in the DB",true);
    op.addOption("structureDB","The path to a structureDB file. Required to verify that the matches are superimposable with the query",true);
    op.addOption("dCut","Cutoff applied when comparing CA distances (default: 0.5 Å");
    op.addOption("RMSDCut","Cutoff applied when calculating RMSD between superimposed termini (default 0.5 Å");
    op.addOption("resectLength","When searching seed termini, allow removal of terminal seed residues up to specified number. This allows for identification of bridges between non-terminal residues. (default = 0)",false);
    op.addOption("minSeedLength","Seeds must retain at least this many residues (only applies when resectLength is > 0)");
    op.addOption("seedA","The path to a seed structure. If provided, will search for overlaps between this seed and other seeds OR seedB only, if that argument is provided",false);
    op.addOption("seedB","The path to a seed structure. If provided, will search for overlaps between this seed and other seeds OR seedA only, if that argument is provided",false);
    op.addOption("seedAPDBList","The path to a list of seed PDB files",false);
    op.addOption("seedBPDBList","The path to a list of seed PDB files",false);
    op.addOption("avoidClashesToStructure","Path to PDB file. If provided, will check for clash between briges and PDB and exclude bridges that clash",false);
    op.addOption("ignoreSeedClash","",false);
    op.addOption("maxSeedTerminiDistance","(default: 15 Å)",false);
    op.addOption("workerID","The 0-indexed unique index for this job, e.g. ${SLURM_ARRAY_TASK_ID} (default = 0)",false);
    op.addOption("nWorkers","The total number of jobs running in parallel (default = 1)",false);
    op.addOption("verbose","");
    // op.addOption("selectedSeed","Either the path to PDB file or the name of a seed in the binary file. If provided, will ONLY search for bridges between this seed and the other seeds (either all seeds in the binary file, or just those in the list)",false);
    op.setOptions(argc,argv);

    string seedBinPath = op.getString("seedBin","");
    string seedListPath = op.getString("seedList","");
    string seedPDBListPath = op.getString("seedPDBList","");
    string bridgeDBPath = op.getString("bridgeDB");
    string structureDBPath = op.getString("structureDB");
    mstreal dCut = op.getReal("dCut",0.5);
    mstreal RMSDCut = op.getReal("RMSDCut",0.5);
    int resectLength = op.getInt("resectLength",0);
    string seedAPath = op.getString("seedA","");
    string seedBPath = op.getString("seedB","");
    string seedAPDBListPath = op.getString("seedAPDBList","");
    string seedBPDBListPath = op.getString("seedBPDBList","");
    string clashStructure = op.getString("avoidClashesToStructure","");
    bool ignoreSeedClash = op.isGiven("ignoreSeedClash");
    int maxSeedTerminiDistance = op.getInt("maxSeedTerminiDistance",15);
    int workerID = op.getInt("workerID",0);
    int nWorkers = op.getInt("nWorkers",1);
    bool verbose = true; // op.isGiven("verbose");

    int overlapLength = 3;
    int minSeedLength = op.getInt("minSeedLength",overlapLength);

    // if ((minimumSeedLength != 0)&(minimumSeedLength < overlapLength)) MstUtils::error("Minimum seed length :"+MstUtils::toString(minimumSeedLength)+" must be at least as great as overlap length: "+MstUtils::toString(overlapLength));
    if ((seedListPath != "")&(seedBinPath=="")) MstUtils::error("If a seed list provided, a seed binary file must also be provided");

    // // Load the seeds from the binary file
    // seedBinaryFile* seeds = nullptr;
    // if (seedBinPath != "") {
    //     seeds = new seedBinaryFile(seedBinPath);
    //     seeds->scanFilePositions();
    // }
    
    // Load the bridge distances/structures from which the distances were calculated
    findSeedBridge bridgeFinder(op.getString("bridgeDB"),MstUtils::toString(workerID));
    bridgeFinder.setMaxSeedDistance(maxSeedTerminiDistance);
    bridgeFinder.loadStructureDB(op.getString("structureDB"));
    bridgeFinder.reportAPVBoundaries();
    
    seedPairDistributor seedPairs(workerID,nWorkers,seedBinPath);
    // By default, perform all-to-all comparison
    if (seedListPath != "") {
        vector<string> seedList = MstUtils::fileToArray(seedListPath);
        seedPairs.setSeeds(seedList);
    }
    if (seedPDBListPath != "") {
        vector<string> seedPDBPathList = MstUtils::fileToArray(seedPDBListPath);
        vector<Structure*> seeds;
        for (string path : seedPDBPathList) {
            Structure* seed = new Structure(path);
            seed->setName(seedPathToName(seed->getName()));
            seeds.push_back(seed);
        }
        seedPairs.setSeeds(seeds);
    }
    if (seedAPDBListPath != "") {
        vector<string> seedPDBPathList = MstUtils::fileToArray(seedAPDBListPath);
        vector<Structure*> seeds;
        for (string path : seedPDBPathList) {
            Structure* seed = new Structure(path);
            seed->setName(seedPathToName(seed->getName()));
            seeds.push_back(seed);
        }
        seedPairs.setSeedsGroupA(seeds);
    }
    if (seedBPDBListPath != "") {
        vector<string> seedPDBPathList = MstUtils::fileToArray(seedBPDBListPath);
        vector<Structure*> seeds;
        for (string path : seedPDBPathList) {
            Structure* seed = new Structure(path);
            seed->setName(seedPathToName(seed->getName()));
            seeds.push_back(seed);
        }
        seedPairs.setSeedsGroupB(seeds);
    }
    if (seedAPath != "") {
        Structure* seedA = new Structure(seedAPath);
        seedA->setName(seedPathToName(seedA->getName()));
        seedPairs.setSeedsGroupA({seedA});
    }
    if (seedBPath != "") {
        Structure* seedB = new Structure(seedBPath);
        seedB->setName(seedPathToName(seedB->getName()));
        seedPairs.setSeedsGroupB({seedB});
    }


    // string bridgeDirName = "seedBridgeStructures/";
    // if (!MstSys::isDir(bridgeDirName)) MstSys::cmkdir(bridgeDirName);
    // string fusedDirName = "fusedSeedBridgeStructures/";
    // if (!MstSys::isDir(fusedDirName)) MstSys::cmkdir(fusedDirName);

    MstTimer timer; timer.start();
    int distanceMatches = 0;
    int rmsdMatches = 0;
    string filename = "";
    string seed_names = "";
    // int minSeedLengthA = minimumSeedLength;
    // int minSeedLengthB = minimumSeedLength;

    fuseSeedsAndBridge fuser(bridgeFinder.bridgeData.getTerminusLength(),MstUtils::toString(workerID));
    if (clashStructure != "") {
        Structure S(clashStructure);
        fuser.setClashCheck(S);
    }

    clashChecker clashCheck;
    cout << "Searching for connections between as many as " << seedPairs.batchSize() << " pairs of seeds" << endl;
    int count = 0;
    while (seedPairs.hasNext()) {
        pair<Structure*,Structure*> seed_pair = seedPairs.next();
        Structure* seedA = seed_pair.first;
        Structure* seedB = seed_pair.second;
        cout << count++ << " " << seedA->getName() << " and " << seedB->getName() << endl;
        if (seedA == seedB) continue;
        // Varying the offset allows us to "resect" the termini in case non-terminal residues are better for forming the bridge
        for (int seedAOffset = 0; seedAOffset <= min(resectLength,max(0,seedA->residueSize()-minSeedLength)); seedAOffset++) {
            if (verbose) cout << "seed A offset: " << seedAOffset << endl;
            vector<Residue*> resectedSeedARes = getSeedResWithResection(seedA,false,seedAOffset);
            if (!ignoreSeedClash) clashCheck.setResidues(resectedSeedARes);
            for (int seedBOffset = 0; seedBOffset <= min(resectLength,max(0,seedB->residueSize()-minSeedLength)); seedBOffset++) {
                if (verbose) cout << "seed B offset: " << seedBOffset << endl;
                
                // Check for clash between seeds (with consideration of the offset)
                if (!ignoreSeedClash) {
                    vector<Residue*> resectedSeedBRes = getSeedResWithResection(seedB,true,seedBOffset);
                    if (clashCheck.checkForClashesToStructure(resectedSeedBRes)) {
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
                fuser.setSeeds(seedA,seedB,seedAOffset,seedBOffset);
                fuser.setBridgeStructures(bridgeFinder.getVerifiedBridgeStructures());
                fuser.writeFusedStructuresToPDB();
            }
        }
    }
    cout << "Done!" << endl;
}