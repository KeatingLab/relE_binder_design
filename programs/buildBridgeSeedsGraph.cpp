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
        return ((group_a_current_index<group_a_max_index)&(group_b_current_index<group_b_max_index));
    }

    pair<Structure*,Structure*> next() {
        if ((group_a_current_index == -1)|(group_b_current_index == -1)) {
            // we round up when determining the batch size, so the final batch is smaller
            int batch_size = ceil(mstreal(seed_group_a.size() * seed_group_b.size()) / nWorkers);
            group_a_current_index = int(workerID * batch_size) / seed_group_a.size();
            group_b_current_index = int(workerID * batch_size) % seed_group_b.size();
            group_a_min_index = group_a_current_index;
            group_b_min_index = group_b_current_index;
            group_a_max_index = min(int(((workerID + 1) * batch_size) / seed_group_a.size()),int(seed_group_a.size()));
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

int main (int argc, char *argv[]) {
    MstOptions op;
    op.setTitle("Given a list of seeds, reports all pairs of seeds that can be joined by a bridging-fragment.");
    op.addOption("seedBin","The path to a seed binary file",false);
    op.addOption("seedList","The path to the list of seeds that will be considered when searching for pairs of connectable seeds",false);
    op.addOption("bridgeDB","The path to a bridgeDB file. If provided, will load and search for matches to the terminus in the DB",true);
    op.addOption("structureDB","The path to a structureDB file. Required to verify that the matches are superimposable with the query",true);
    op.addOption("dCut","Cutoff applied when comparing CA distances (default: 0.5 Å");
    op.addOption("RMSDCut","Cutoff applied when calculating RMSD between superimposed termini (default 0.5 Å");
    op.addOption("resectLength","When searching seed termini, allow removal of terminal seed residues up to specified number. This allows for identification of bridges between non-terminal residues. (default = 0)",false);
    op.addOption("minSeedLength","Seeds must retain at least this many residues (only applies when resectLength is > 0)");
    op.addOption("seedA","The path to a seed structure. If provided, will search for overlaps between this seed and other seeds OR seedB only, if that argument is provided",false);
    op.addOption("seedB","The path to a seed structure. If provided, will search for overlaps between this seed and other seeds OR seedA only, if that argument is provided",false);
    op.addOption("avoidClashesToStructure","Path to PDB file. If provided, will check for clash between briges and PDB and exclude bridges that clash",false);
    // op.addOption("selectedSeed","Either the path to PDB file or the name of a seed in the binary file. If provided, will ONLY search for bridges between this seed and the other seeds (either all seeds in the binary file, or just those in the list)",false);
    op.setOptions(argc,argv);

    string seedBinPath = op.getString("seedBin","");
    string seedListPath = op.getString("seedList","");
    string bridgeDBPath = op.getString("bridgeDB");
    string structureDBPath = op.getString("structureDB");
    mstreal dCut = op.getReal("dCut",0.5);
    mstreal RMSDCut = op.getReal("RMSDCut",0.5);
    int resectLength = op.getInt("resectLength",0);
    string seedAPath = op.getString("seedA","");
    string seedBPath = op.getString("seedB","");
    string clashStructure = op.getString("avoidClashesToStructure","");

    int overlapLength = 3;
    int minSeedLength = op.getInt("minSeedLength",overlapLength);

    // if ((minimumSeedLength != 0)&(minimumSeedLength < overlapLength)) MstUtils::error("Minimum seed length :"+MstUtils::toString(minimumSeedLength)+" must be at least as great as overlap length: "+MstUtils::toString(overlapLength));
    if ((seedListPath != "")&(seedBinPath=="")) MstUtils::error("If a seed list provided, a seed binary file must also be provided");

    // Load the seeds from the binary file
    seedBinaryFile* seeds = nullptr;
    if (seedBinPath != "") {
        seeds = new seedBinaryFile(seedBinPath);
        seeds->scanFilePositions();
    }
    
    // Load the bridge distances/structures from which the distances were calculated
    findSeedBridge bridgeFinder(op.getString("bridgeDB"));
    bridgeFinder.loadStructureDB(op.getString("structureDB"));
    bridgeFinder.reportAPVBoundaries();

    // Load the seeds which we will search for bridges between
    vector<Structure*> all_loaded_seeds;
    Structure* selected_seedA = nullptr;
    Structure* selected_seedB = nullptr;
    vector<Structure*> selectedSeedsGroupA;
    vector<Structure*> selectedSeedsGroupB;
    if (seedAPath != "") {
        selected_seedA = new Structure(seedAPath);
        selected_seedA->setName(MstSys::splitPath(selected_seedA->getName(),1));
        all_loaded_seeds.push_back(selected_seedA);
    }
    if (seedBPath != "") {
        selected_seedB = new Structure(seedBPath);
        selected_seedB->setName(MstSys::splitPath(selected_seedB->getName(),1));
        all_loaded_seeds.push_back(selected_seedB);
    }
    // By default, perform all-to-all comparison
    if (seedListPath != "") {
        vector<string> seedList = MstUtils::fileToArray(seedListPath);
        for (string seedName : seedList) {
            Structure* seed = seeds->getStructureNamed(seedName);
            selectedSeedsGroupA.push_back(seed);
            selectedSeedsGroupB.push_back(seed);
            all_loaded_seeds.push_back(seed);
        }
    }
    // If seed A and/or seed B were specified, then perform one-to-one/one-to-all comparison
    if (selected_seedA != nullptr) selectedSeedsGroupA = {selected_seedA};
    if (selected_seedB != nullptr) selectedSeedsGroupB = {selected_seedB};

    string bridgeDirName = "seedBridgeStructures/";
    if (!MstSys::isDir(bridgeDirName)) MstSys::cmkdir(bridgeDirName);
    // string fusedDirName = "fusedSeedBridgeStructures/";
    // if (!MstSys::isDir(fusedDirName)) MstSys::cmkdir(fusedDirName);

    MstTimer timer; timer.start();
    int distanceMatches = 0;
    int rmsdMatches = 0;
    string filename = "";
    string seed_names = "";
    // int minSeedLengthA = minimumSeedLength;
    // int minSeedLengthB = minimumSeedLength;

    fuseSeedsAndBridge fuser(bridgeFinder.bridgeData.getTerminusLength());
    if (clashStructure != "") {
        Structure S(clashStructure);
        fuser.setClashCheck(S);
    }
    clashChecker clashCheck;
    cout << "Searching for connections between " << selectedSeedsGroupA.size() << " x " << selectedSeedsGroupB.size() << " seeds" << endl;
    for (int i = 0; i < selectedSeedsGroupA.size(); i++) {
        Structure* seedA = selectedSeedsGroupA[i];
        clashCheck.setStructure(*seedA);
        for (int j = 0; j < selectedSeedsGroupB.size(); j++) {
            Structure* seedB = selectedSeedsGroupB[j];
            if (seedA == seedB) continue;
            if ((selectedSeedsGroupA.size()>1)&&(selectedSeedsGroupB.size()>1)&&(j > i)) continue;
            if (clashCheck.checkForClashesToStructure(seedB->getResidues())) continue;
            // Varying the offset allows us to "resect" the termini in case non-terminal residues are better for forming the bridge
            for (int seedAOffset = 0; seedAOffset <= min(resectLength,max(0,seedA->residueSize()-minSeedLength)); seedAOffset++) {
                for (int seedBOffset = 0; seedBOffset <= min(resectLength,max(0,seedB->residueSize()-minSeedLength)); seedBOffset++) {
                    // seed A -> seed B
                    bridgeFinder.setSearchQuery(seedA,seedB,seedAOffset,seedBOffset);
                    timer.start();
                    distanceMatches = bridgeFinder.searchSeedsByCADistance(dCut);
                    timer.stop();
                    std::cout << "Took " << timer.getDuration(MstTimer::msec) << " ms to find matches" << endl;
                    if (distanceMatches == 0) continue;
                    timer.start();
                    rmsdMatches = bridgeFinder.verifyMatchesBySuperposition(RMSDCut);
                    timer.stop();
                    std::cout << "Took " << timer.getDuration(MstTimer::msec) << " ms to verify the matches" << endl;
                    if (rmsdMatches == 0) continue;
                    cout << "Found bridges between " << seedA->getName() << " with offset " << MstUtils::toString(seedAOffset) << " and " << seedB->getName() << " with offset " << MstUtils::toString(seedBOffset) << endl;

                    // Save bridge structure info
                    seed_names = seedA->getName()+"_"+MstUtils::toString(seedAOffset)+"_"+MstUtils::toString(seedBOffset)+"_"+seedB->getName();
                    bridgeFinder.getVerifiedBridgeLengthDist(seed_names);
                    bridgeFinder.writeToMultiPDB(bridgeDirName+"/"+seed_names+".pdb",seed_names);

                    // Fuse seeds with bridge
                    fuser.setSeeds(seedA,seedB,seedAOffset,seedBOffset);
                    fuser.setBridgeStructures(bridgeFinder.getVerifiedBridgeStructures());
                    fuser.writeFusedStructuresToPDB();

                    // seed B -> seed A
                    bridgeFinder.setSearchQuery(seedB,seedA,seedBOffset,seedAOffset);
                    timer.start();
                    distanceMatches = bridgeFinder.searchSeedsByCADistance(dCut);
                    timer.stop();
                    std::cout << "Took " << timer.getDuration(MstTimer::msec) << " ms to find matches" << endl;
                    if (distanceMatches == 0) continue;
                    timer.start();
                    rmsdMatches = bridgeFinder.verifyMatchesBySuperposition(RMSDCut);
                    timer.stop();
                    std::cout << "Took " << timer.getDuration(MstTimer::msec) << " ms to verify the matches" << endl; 
                    if (rmsdMatches == 0) continue;
                    cout << "Found bridges between " << seedB->getName() << " with offset " << MstUtils::toString(seedBOffset) << " and " << seedA->getName() << " with offset " << MstUtils::toString(seedAOffset) << endl;

                    // Save bridge structure info
                    seed_names = seedB->getName()+"_"+MstUtils::toString(seedBOffset)+"_"+MstUtils::toString(seedAOffset)+"_"+seedA->getName();
                    bridgeFinder.getVerifiedBridgeLengthDist(seed_names);
                    bridgeFinder.writeToMultiPDB(bridgeDirName+"/"+seed_names+".pdb",seed_names);

                    // Fuse seeds with bridge
                    fuser.setSeeds(seedB,seedA,seedBOffset,seedAOffset);
                    fuser.setBridgeStructures(bridgeFinder.getVerifiedBridgeStructures());
                    fuser.writeFusedStructuresToPDB();
                }
            }
        }
    }

    if (seeds != nullptr) delete seeds;
    for (Structure* seed : all_loaded_seeds) delete seed;

    cout << "Done!" << endl;
}