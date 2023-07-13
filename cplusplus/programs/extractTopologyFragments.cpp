#include "mstfuser.h"
#include "mstrotlib.h"
#include "msttypes.h"
#include "mstoptions.h"

#include "bridgeseeds.h"
#include "generateseeds.h"
#include "residueframe.h"
#include "utilities.h"
#include "utilitiesio.h"

int main (int argc, char *argv[]) {
    MstOptions op;
    op.setTitle("Load topo, save fragments to multiPDB");
    op.addOption("previousRunDir","The path to a directory where connectFragments was previously run. Provide this if you're connecting fragments that have already been fused and you want to preserve their topology information",false);
    op.addOption("topoName","Name of the topology to extract",true);
    // op.addOption("selectedSeed","Either the path to PDB file or the name of a seed in the binary file. If provided, will ONLY search for bridges between this seed and the other seeds (either all seeds in the binary file, or just those in the list)",false);
    op.setOptions(argc,argv);

    MstTimer timer; timer.start();

    string previousRunDirPath = op.getString("previousRunDir",".");
    string topoName = op.getString("topoName");

    // load old topoDB
    string oldTopoDBpath = previousRunDirPath + "/topoDB.json";
    topologyDB oldTopoDB;
    oldTopoDB.readDBFromJSON(oldTopoDBpath);

    // load old fragDB
    string oldFragDBpath = previousRunDirPath + "/fragmentDB.pdb";
    multiPDBFile oldFragDB(oldFragDBpath);
    oldFragDB.countStructures();

    if (!oldTopoDB.isTopoInDB(topoName)) MstUtils::error("topology not in DB");

    shared_ptr<fragmentTopology> fT = oldTopoDB.getTopologyFromDB(topoName,oldFragDB);

    multiPDBFile newFragDB(topoName+"_fragmentDB.pdb",false);
    fT->addFragmentsInTopoToDB(newFragDB);

    cout << "Done!" << endl;
}