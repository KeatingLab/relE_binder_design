#include "mstsystem.h"
#include "mstoptions.h"

#include "residuecontact.h"

int main(int argc, char *argv[]) {
    MstOptions op;
    op.setTitle("For testing the vdwContacts method.");
    op.addOption("p", "The structure which will be used to compute contact information",true);
    op.setOptions(argc,argv);
    
    Structure S(op.getString("p"));
    string structure_name = MstSys::splitPath(S.getName(),1);
    vdwContacts vdwC(&S);
    
    fstream out;
    MstUtils::openFile(out, structure_name+"_vdwContacts.csv",fstream::out);
    out << "structure,resi_chainid,resi_num,resi_name,resj_chainid,resj_num,resj_name" << endl;

    fstream cont_out;
    MstUtils::openFile(cont_out, structure_name+"_drawcontacts.tsv",fstream::out);
    
    for (Residue* Ri : S.getResidues()) {
//        cout << "res: " << R->getNum() << endl;
        set<Residue*> contactingResidues = vdwC.getInteractingRes(Ri);
        for (Residue* Rj : contactingResidues) {
            out << S.getName() << ",";
            out << Ri->getChainID() << "," << Ri->getNum() << "," << Ri->getName() << ",";
            out << Rj->getChainID() << "," << Rj->getNum() << "," << Rj->getName() << ",";
            out << endl;

            cont_out << Ri->getChainID() << Ri->getNum() << "\t";
            cont_out << Rj->getChainID() << Rj->getNum() << "\t";
            cont_out << 1;
            cont_out << endl;
        }
    }
    out.close();
    cont_out.close();

    return 0;
}
