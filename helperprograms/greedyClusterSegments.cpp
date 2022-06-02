#include "mstrotlib.h"
#include "msttypes.h"
#include "mstoptions.h"
#include "mstfasst.h"

// #include "alignframes.h"
// #include "residuecontact.h"
// #include "residueframe.h"

#include "fragmentdb.h"

int main (int argc, char *argv[]) {
    MstOptions op;
    op.setTitle("Loads PDB dataset, greedily clusters segments, and writes a new database file");
    op.addOption("fasstDB","A path to a structure database that is compatible with FASST. Note: DB should have been created with '--c', so that there are no gaps within a chain",true);
    op.addOption("segLength","The length (residues) of the segments extracted from structures (default: 8)");
    op.addOption("clusterRMSDCutoff","The RMSD cutoff (Å) that is used to define greedy clusters (default:  0.25)");
    op.addOption("coverage","The fraction of segments that should be covered when clustering (default 0.9)");
    op.addOption("nSubsampled","The number of segments to sample when selecting a cluster centroid (default 1000)");
    op.setOptions(argc,argv);

    string fasstDBPath = op.getString("fasstDB");
    int segLength = op.getInt("segLength",4);
    mstreal clusterRMSDCutoff = op.getReal("clusterRMSDCutoff",0.5);
    mstreal coverage = op.getReal("coverage",0.9);
    int nSubsampled = op.getInt("nSubsampled",1000);

    clusterDBSegments clusterer(segLength);
    
    clusterer.loadFASSTDB(fasstDBPath);

    clusterer.cluster(clusterRMSDCutoff,coverage,nSubsampled);

    clusterer.writeClusterRepresentativesToBinFile("clusterRepresentatives.bin");
    clusterer.writeClusterRepresentativesToPDBFile("top100clusterRepresentatives.pdb",100);
    clusterer.writeClusterMembersToPDBFile("top10Clusters.pdb",10,10);
    clusterer.writeClusterInfo("clusterInfo.csv");

    cout << "Done!" << endl;
}