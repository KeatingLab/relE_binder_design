#include "mstrotlib.h"
#include "msttypes.h"
#include "mstoptions.h"
#include "mstfasst.h"

#include "fragmentdb.h"
// #include "alignframes.h"
// #include "residuecontact.h"
// #include "residueframe.h"

int main (int argc, char *argv[]) {
    MstOptions op;
    op.setTitle("Calculates the designability (number of matches) for library of segments and then finds overlaps between the termini of protein backbone segments");
    op.addOption("segmentBin","The path to a binary file containing protein backbone segments",true);
    op.addOption("fasstDB","The path to a FASST DB that is search to estimate the designability of each segment. If not provided, will skip the search",false);
    op.addOption("segmentMatchRMSD","If FASST DB is provided, can optionally specify the RMSD used to search for matches to each segment in the DB (default 0.5 Å)",false);
    op.addOption("overlapRMSD","The RMSD cutoff (Å) that is used to define overlap between the segments (default:  0.5)");
    op.addOption("name","The name of the segment overlap graph",true);
    op.setOptions(argc,argv);

    string segmentBinPath = op.getString("segmentBin");
    string fasstDBPath = op.getString("fasstDB","");
    mstreal segmentMatchRMSD = op.getReal("segmentMatchRMSD",0.5);
    mstreal overlapRMSD = op.getReal("overlapRMSD",0.5);
    string name = op.getString("name");

    segmentGraph graph;

    // Load segments from binary file
    vector<Structure*> segments = clusterDBSegments::readClusterRepresentativesFromBinFile(segmentBinPath);

    // If fasst DB is provided, search segments against the database to get the number of matches
    vector<int> segmentsMatchNumber;
    if (fasstDBPath != "") {
        FASST F;
        F.readDatabase(fasstDBPath,2); // Only need number of matches, so should be able to load in mode 2
        fasstSearchOptions op;
        op.setRMSDCutoff(segmentMatchRMSD);
        cout << "Searching segments for matches" << endl;
        for (Structure* segment: segments) {
            cout << "Segment: " << segment->getName() << endl;
            F.setQuery(*segment);
            F.search();
            int numMatches = F.numMatches();
            cout << "Found " << numMatches << endl;
            segmentsMatchNumber.push_back(numMatches);
        }
    }

    // Construct graph
    if (fasstDBPath == "") graph.constructGraphFromSegments(segments,overlapRMSD);
    else graph.constructGraphFromSegments(segments,overlapRMSD,segmentsMatchNumber);

    graph.writeGraphFile(name);

    cout << "Done!" << endl;
}