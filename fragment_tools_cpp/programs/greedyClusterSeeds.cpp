#include "msttypes.h"
#include "mstoptions.h"

#include "alignframes.h"
#include "generateseeds.h"
#include "residuecontact.h"
#include "residueframe.h"
#include "utilitiesio.h"

#include <chrono>

int main(int argc, char *argv[])
{
    MstOptions op;
    op.setTitle("Greedily clusters seeds.");
    op.addOption("binderBin","Path to the seed binary file",true);
    op.addOption("seedList","Path to a file where each line is the name of a seed in the binary file. If provided, will only cluster these seeds");
    op.addOption("clusterRMSDCutoff","The RMSD cutoff (Å) that is used to define greedy clusters (default:  0.5)");
    op.addOption("seedLength","The number of residues in each seed (default 5)");
    op.addOption("coverage","Greedy clustering will continue until this fraction of seeds are covered. The remaining seeds will be placed in their own clusters (default 0.5)");
    op.setOptions(argc,argv);

    string binderBinPath = op.getString("binderBin");
    string seedListPath = op.getString("seedList","");
    mstreal clusterRMSDCutoff = op.getReal("clusterRMSDCutoff",0.5);
    int seedLength = op.getInt("seedLength",5);
    mstreal coverage = op.getReal("coverage",0.5);

    set<string> selectedSeeds;
    if (seedListPath != "") {
        vector<string> seedList = MstUtils::fileToArray(seedListPath);
        selectedSeeds = set<string>(seedList.begin(),seedList.end());
    }

    /*
    Most seeds will be the same length, but some will be shorter. Gather all seeds of the max length, and
    cluster those. Then assign the remaining seeds to the clusters, or make new ones.
    
    jk for now throw out everything shorter than five residues
    */

    cout << "Extracting seeds from binary file" << endl;
    vector<Structure*> seeds;
    vector<vector<Atom*>> seedAtoms;
    seedBinaryFile seedBin(binderBinPath);
    // iterate over file and get seeds with correct length that are in the lists
    while (seedBin.hasNext()) {
        Structure* seed = seedBin.next();
        if (selectedSeeds.count(seed->getName()) == 0) {
            delete seed;
            continue;
        }
        if (seed->residueSize() != seedLength) {
            delete seed;
            continue;
        }
        seeds.push_back(seed);
        seedAtoms.push_back(seed->getAtoms());
    }
    cout << "Extracted " << seeds.size() << " seeds total from the file" << endl;

    int Nmax = 1000;
    Clusterer Cl;
    Cl.optimizeAlignments(false);
    cout << "Greedily clustering the segments at an RMSD cutoff of " << clusterRMSDCutoff << endl;
    vector<vector<int>> clusters = Cl.greedyCluster(seedAtoms,clusterRMSDCutoff,Nmax,coverage,-1,true);
    cout << "Defined " << clusters.size() << " clusters" << endl;

    cout << "Writing cluster info" << endl;
    fstream info_out;
    MstUtils::openFile(info_out,"clusterInfo.csv",fstream::out);
    info_out << "cluster_ID,seedName,representative" << endl;
    set<int> clusteredSeedIdx;
    bool rep = false;
    for (int i = 0; i < clusters.size(); i++) {
        for (int j = 0; j < clusters[i].size(); j++) {
            rep = false;
            if (j == 0) rep = true;
            info_out << i << "," << seeds[clusters[i][j]]->getName() << ",";
            info_out << rep << endl;
            clusteredSeedIdx.insert(clusters[i][j]);
        }
    }
    int clusterIdx = clusters.size();
    for (int i = 0; i < seeds.size(); i++) {
        if (clusteredSeedIdx.count(i) > 0) continue;
        info_out << clusterIdx << "," << seeds[i]->getName() << ",1" << endl;
        clusterIdx++;
    }
    info_out.close();

    // cout << "Writing PDB files of representatives for top 100 clusters" << endl;
    // fstream pdb_out;
    // MstUtils::openFile(pdb_out,"clusterRepresentativestop100.pdb",fstream::out);
    // int count = 0;
    // for (int i = 0; i < clusters.size(); i++) {
    //     if (i >= 100) break;
    //     string name = "cluster_"+MstUtils::toString(i);
    //     pdb_out << "HEADER     " << name << endl;
    //     Structure* seed = seeds[clusters[i][0]];
    //     seed->writePDB(pdb_out);

    //     count += clusters[i].size();
    // }
    // pdb_out.close();

    return 0;
}