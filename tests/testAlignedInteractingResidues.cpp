#include "mstsystem.h"
#include "mstoptions.h"
#include "mstsequence.h"

#include "alignframes.h"
#include "residuecontact.h"
#include "residueframe.h"

int main(int argc, char *argv[]) {
    MstOptions op;
    op.setTitle("For testing the aligned interacting residues method");
    op.addOption("p", "The structure from which interacting residues will be extracted",true);
    op.setOptions(argc,argv);
    
    Structure S(op.getString("p"));
    string structure_name = MstSys::splitPath(S.getName(),1);
    vdwContacts vdwC(S.getResidues());

    alignInteractingFrames alignF;
    proteinFrameDB& db = alignF.getDB();
    augmentedStructure aS(S);
    db.addTarget(aS);
    
    db.setVDWContacts(0,vdwC.getAllInteractingRes());

    residueFrame rF;
    alignF.setRefFrame(rF);

    set<string> all_residue_aa;
    for (Residue* R : aS.getResidues()) all_residue_aa.insert(R->getName());

    cout << "Structure has " << all_residue_aa.size() << " amino acid types" << endl;
    for (string residue_aa : all_residue_aa) {
        cout << "aa: " << residue_aa << endl;

        alignF.setAA(residue_aa);
        alignF.findMobileFrames();
        cout << "Structure has " << alignF.getNumInteracting() << " residue interactions from this residue type" << endl;

        alignF.writeAlignedInteractingResToPDB(residue_aa+"_interactions.pdb");
    }
    return 0;
}