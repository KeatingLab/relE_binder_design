#include "mstsequence.h"
#include "mstsystem.h"
#include "mstoptions.h"

#include "alignframes.h"
#include "hashframes.h"
#include "residuecontact.h"
#include "residueframe.h"

void printFrameInfo(residueFrame* frame) {
    anglesFromFrame angles(frame);
    cout << "x: " << frame->getO()[0];
    cout << "\ty: " << frame->getO()[1];
    cout << "\tz: " << frame->getO()[2];
    cout << "\talpha: " << 360*angles.getAlpha()/(2*M_PI);
    cout << "\tbeta: " << 360*angles.getBeta()/(2*M_PI);
    cout << "\tgamma: " << 360*angles.getGamma()/(2*M_PI);
    cout << endl;
}

int main(int argc, char *argv[]) {
    MstOptions op;
    op.addOption("pdb","the structure from which a residue pair will be scored");
    op.addOption("sel","a selection consisting of two residues in the provided structure");
    op.addOption("frame_db","Path to db");
    op.addOption("struct_db","Path to db of structures");
    op.setOptions(argc,argv);

    MstTimer timer;

    fstream out;
    MstUtils::openFile(out,"residueFrames.tsv",fstream::out);

    // extract the residue pair query from the structure
    augmentedStructure S(op.getString("pdb"));
    selector sel(S);
    vector<Residue*> query_res = sel.selectRes(op.getString("sel"));
    MstUtils::assert(query_res.size(),"The selection must consist of two residues only","testSearchHashTable::main");
    cout << "Selected residues: " << query_res[0]->getChainID() << query_res[0]->getNum() << " and " << query_res[1]->getChainID() << query_res[1]->getNum() << endl;
    residueFrame* frameI = S.getResidueFrame(query_res[0]->getResidueIndex());
    residueFrame* frameJ = S.getResidueFrame(query_res[1]->getResidueIndex());
    residueFrame* query = frameJ->frameRelativeToOther(*frameI);

    // // find the transformation between other and the global reference frame
    // Transform other_to_ref = TransformFactory::switchFrames(Frame(), *frameI);

    // // apply the transformation to this frame and return
    // residueFrame* query = new residueFrame(*frameJ);
    // other_to_ref.apply(*query);
    // Structure transformedStructure(S);
    // other_to_ref.apply(transformedStructure);
    // transformedStructure.writePDB("transformedStructure.pdb");
    
    residueFrame ref;

    frameJ->writeToFile("queryFrameJ",out);
    frameI->writeToFile("queryFrameI",out);
    query->writeToFile("queryFrame_inglobalframe",out);
    ref.writeToFile("globalRefFrame",out);

    //write original residue pair
    Structure original_res_pair(query_res);
    original_res_pair.writePDB("original_res_pair.pdb");
    
    //write residue pair after transformation
    // other_to_ref.apply(original_res_pair);
    original_res_pair.writePDB("original_res_pair_transformed.pdb");

    res_t resI_aa = SeqTools::aaToIdx(query_res[0]->getName());

    // load the frames from the DB into the hash table
    string dbPath = op.getString("frame_db");
    frameDB DB(dbPath);

    vector<mobileFrame*> loaded;
    while (DB.hasNext()) {
        mobileFrame* mF = DB.next();

        if (mF->getResIIndex() == resI_aa) loaded.push_back(mF);
        else delete mF;
    }
    cout << "Loaded " << loaded.size() << " frames interacting with residue type: " << query_res[0]->getName() << endl;

    // get the matches to the query
    boundingBox bbox;
    bool verbose = true;
    bbox.update(loaded);
    frameTable table(bbox, 0.25, 36);

    for (mobileFrame* frame : loaded) {
        table.insertFrame(frame);
    }
    cout << "Done loading frames into table" << endl;
    timer.start();
    set<mobileFrame*> matches = table.findSimilarFrames(query, 1.0, 45.0);
    timer.stop();
    cout << "Found " << matches.size() << " frames matching the criteria in " << timer.getDuration(MstTimer::timeUnits::msec) << " ms" << endl;

    cout << "Query frame: " << endl;
    printFrameInfo(query);

    // extract the interacting residues from the DB
    alignInteractingFrames alignF(op.getString("struct_db"));
    Transform tf; int count = 0;
    for (mobileFrame* frame : matches) {
        frame->writeToFile("matchFrame_"+MstUtils::toString(count),out);
        cout << "Match frame: " << endl;
        printFrameInfo(frame);

        Structure* res = alignF.getAlignedInteractingRes(frame);
        res->writePDB(res->getName()+"_inglobalframe.pdb");
        // align to the original query
        tf = TransformFactory::switchFrames(*frameI,Frame());
        tf.apply(res);
        res->writePDB(res->getName()+".pdb");
        delete res;
        count++;
    }
    delete query;

    return 0;
}