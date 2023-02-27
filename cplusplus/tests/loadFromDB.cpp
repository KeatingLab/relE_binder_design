#include "msttypes.h"
#include "mstoptions.h"

#include "alignframes.h"
#include "residuecontact.h"
#include "residueframe.h"
#include "utilities.h"

#include <chrono>

int main(int argc, char *argv[]) {
    MstOptions op;
    op.setTitle("Creates a mobile frame database from a protein structure database.");
    op.addOption("db","Path to protein structure database");
    op.addOption("target_id","");
    op.setOptions(argc,argv);

    augmentedStructureDB DB(op.getString("db"));
    const augmentedStructure& S = DB.getTarget(op.getInt("target_id"));
    string pdb_path = "target_"+MstUtils::toString(op.getInt("target_id"))+".pdb";
    S.writePDB(pdb_path);

    cout << "Done!" << endl;
    return 0;
}