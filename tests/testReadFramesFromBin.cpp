#include "mstsystem.h"
#include "mstoptions.h"

#include "alignframes.h"
#include "residuecontact.h"
#include "residueframe.h"

int main(int argc, char *argv[]) {
    MstOptions op;
    op.addOption("struct_db","Path to db");
    op.addOption("frame_db","Path to db");
    op.addOption("rate","subsample rate");
    op.setOptions(argc,argv);

    mstreal rate = op.getReal("rate",0.0001);
    string structDB = op.getString("struct_db");
    alignInteractingFrames alignF(structDB);

    string dbPath = op.getString("frame_db");
    frameDB DB(dbPath);

    while (DB.hasNext()) {
        if (MstUtils::randUnit() > rate) {
            DB.skip();
            continue;
        }
        mobileFrame* mF = DB.next();
        cout << "frame: " << mF->getName() << endl;
        cout << *mF << endl;

        Structure* intRes = alignF.getAlignedInteractingRes(mF);
        intRes->writePDB(mF->getName()+".pdb");

        delete mF;
        delete intRes;
    }
    return 0;
}