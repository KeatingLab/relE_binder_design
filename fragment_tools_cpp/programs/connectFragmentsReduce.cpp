#include "mstfuser.h"
#include "mstrotlib.h"
#include "msttypes.h"
#include "mstoptions.h"

#include "bridgeseeds.h"
#include "generateseeds.h"
#include "residueframe.h"
#include "utilities.h"

int main (int argc, char *argv[]) {
    MstOptions op;
    op.setTitle("After all connectFragmentsMapper workers have completed, run this program to unify the database files");
    op.setOptions(argc,argv);

    // Load info file to see how many jobs were attempted
    vector<string> lines = MstUtils::fileToArray("job_info.txt");
    if (lines.size() != 6) MstUtils::error("Wrong number of lines in info file","main");
    vector<string> split_line = MstUtils::split(lines[5]);
    if (split_line[0] != "nworkers") MstUtils::error("Line is ill-formatted","main");
    int nWorkers = MstUtils::toInt(split_line[1]);

    // Load the DB files (fusedDB, fragmentDB, topoDB) from each worker and add to new, combined file
    multiPDBFile writeFusedDB("fusedDB.pdb",false), writeFragmentDB("fragmentDB.pdb",false);
    topologyDB writeTopoDB;
    string path;
    for (int workerID = 0; workerID < nWorkers; workerID++) {
        cout << "workerID: " << workerID << endl;
        path = "fusedDB_"+MstUtils::toString(workerID)+".pdb";
        multiPDBFile readFusedDB(path);
        cout << readFusedDB.countStructures() <<  " structures in fusedDB" << workerID << endl;
        readFusedDB.copyStructuresToNewMultiPDB(writeFusedDB);

        path = "fragmentDB_"+MstUtils::toString(workerID)+".pdb";
        multiPDBFile readFragmentDB(path);
        cout << readFragmentDB.countStructures() <<  " structures in fragmentDB" << workerID << endl;
        readFragmentDB.copyStructuresToNewMultiPDB(writeFragmentDB);

        path = "topoDB_"+MstUtils::toString(workerID)+".json";
        writeTopoDB.readDBFromJSON(path,true);
    }
    writeTopoDB.writeDBtoJSON("topoDB.json");

    cout << "Done!" << endl;
}