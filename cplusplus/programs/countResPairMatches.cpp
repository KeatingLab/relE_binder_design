#include "msttypes.h"
#include "mstoptions.h"

#include "alignframes.h"
#include "generateseeds.h"
#include "scoreinterface.h"
#include "residuecontact.h"

#include <chrono>

vector<string> getBatch(vector<string> all_paths, int worker, int nWorkers) {
    // Assumes that the jobs are sorted from longest to shortest (e.g. most to fewest residues)
    // Later batches will be smaller if the number of workers is close to the number of batches
    // 1) Pair the longest and shortest jobs (each pair will have close to equal total length)
    int all_jobs = all_paths.size();
    int final_job_idx = 0;
    if (all_jobs % 2 != 0) {
        // odd number, remove the final job add to first batch
        final_job_idx = all_jobs - 1;
        all_jobs--;
    }
    int n_pairs = all_jobs / 2;
    // 2) Assemble batches from pairs
    vector<string> batch_paths;
    int i = worker;
    while (i < n_pairs) {
        if ((i == 0)&&(final_job_idx != 0)) {
            // If it exists, add the "odd one out" to the first batch
            batch_paths.push_back(all_paths[final_job_idx]);
        }
        // first job in pair
        batch_paths.push_back(all_paths[i]);
        // second job in pair
        batch_paths.push_back(all_paths[all_jobs-1-i]);
        i += nWorkers;
    }
    return batch_paths;
}

int main(int argc, char *argv[]) {
    MstOptions op;
    op.setTitle("Defines k-NN around each residue in a protein and searches for matches to all unique residue pairs.");
    op.addOption("resPairDB","Path to residue pair database",true);
    op.addOption("structurePDB","A PDB file defining the structure to be scored",false);
    op.addOption("structurePDBList","A file where each line is the path to a PDB file defining the structure to be scored",false);
    op.addOption("kNeighbors","The number of nearest neighboring residues (including self) to consider when searching for matches (default: 30)",false);
    // op.addOption("RMSDCutoff","The RMSD cutoff that is applied when confirming a putative match (default = 0.25 Å)",false);
    op.addOption("worker","The ID of the worker assigned to this batch (default = 1)");
    op.addOption("nWorkers","The total number of workers (default = 1)");
    op.addOption("noSearch","Skip search and only output kNN residue info");
    op.addOption("verbose","");
    op.setOptions(argc,argv);

    if ((op.isGiven("structurePDB") && op.isGiven("structurePDBList"))||(!op.isGiven("structurePDB") && !op.isGiven("structurePDBList"))) {
        MstUtils::error("Must provide either --structurePDB or --structurePDBList, but not both");
    } else if (op.isGiven("structurePDB")) {
        cout << "structurePDB mode" << endl;
    } else {
        cout << "structurePDBList mode" << endl;
    } 

    string resPairDB = op.getString("resPairDB");
    string structurePDB = op.getString("structurePDB","");
    string structurePDBList = op.getString("structurePDBList","");
    int kNeighbors = op.getInt("kNeighbors",30);
    // mstreal RMSDCutoff = op.getReal("RMSDCutoff",0.25);
    bool verbose = op.isGiven("verbose");
    int worker = op.getInt("worker",1);
    int nWorkers = op.getInt("nWorkers",1);
    bool noSearch = op.getBool("noSearch",false);

    findkNNResPairs kNNresPairs(resPairDB,kNeighbors,noSearch);

    if (structurePDB != "") {
        Structure originalStructure(structurePDB,"SKIPHETERO|ALLOW ILE CD1");
        Structure cleanedStructure;
        RotamerLibrary::extractProtein(cleanedStructure,originalStructure);
        augmentedStructure structure(cleanedStructure);
        structure.setName(MstSys::splitPath(structurePDB,1));
        kNNresPairs.setStructure(&structure);
        kNNresPairs.findkNN();
        kNNresPairs.searchkNN();
        kNNresPairs.writeToFile();

    } else if (structurePDBList != "") {
        vector<string> lines = MstUtils::fileToArray(structurePDBList);
        vector<string> batch;
        if (nWorkers > 1) {
            vector<string> batch = getBatch(lines,worker,nWorkers);

            fstream out;
            MstUtils::openFile(out,MstUtils::toString(worker)+"_batch_list.txt",fstream::out);
            for (string structurePDB : batch) {
                out << structurePDB << endl;
            }
        } else {
            batch = lines;
        }

        for (string structurePDB : batch) {
            string structureName = MstSys::splitPath(structurePDB,1);
            string outputName = structureName+"_respair_matches.csv";
            cout << "Trying to get data for " << structurePDB << endl;
            if (ifstream(outputName)) {
                cout << "Job is done, skipping..." << endl;
                continue;
            }
            Structure originalStructure(structurePDB,"SKIPHETERO|ALLOW ILE CD1");
            Structure cleanedStructure;
            RotamerLibrary::extractProtein(cleanedStructure,originalStructure);
            augmentedStructure structure(cleanedStructure);
            structure.setName(structureName);
            kNNresPairs.setStructure(&structure);
            kNNresPairs.findkNN();
            kNNresPairs.searchkNN();
            kNNresPairs.writeToFile();
        }
    }
    cout << "Done" << endl;

    return 0;
}