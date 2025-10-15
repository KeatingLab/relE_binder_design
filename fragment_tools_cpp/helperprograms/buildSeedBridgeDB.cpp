#include "msttypes.h"
#include "mstoptions.h"

#include "bridgeseeds.h"

int main (int argc, char *argv[]) {
    MstOptions op;
    op.setTitle("Reports geometric parameters that reflect whether two segments can be easily joined by building in intervening residues");
    op.addOption("structureDB","The path to a PDB structure database that will be searched to find seed bridge data",true);
    op.addOption("prefix","The prefix will be added to the seed bridge database file name",false);
    op.addOption("maxLength","The maximum length of bridge to store in the DB (default: 12 residues)");
    op.setOptions(argc,argv);

    string structureDBPath = op.getString("structureDB");
    string prefix = op.getString("prefix","");
    int maxLength = op.getInt("maxLength",12);

    seedBridgeDB bridgeData(maxLength);
    bridgeData.loadProteinStructures(structureDBPath);
    bridgeData.buildDBfromStructures(maxLength);
    bridgeData.writeDBtoFile(prefix+".bridge.bin");

    cout << "Done!" << endl;
}