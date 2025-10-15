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
    vdwContacts vdwC(S.getResidues());
    
    fstream out;
    MstUtils::openFile(out, structure_name+"_vdwContacts.csv",fstream::out);
    out << "structure,resi_chainid,resi_num,resi_name,resj_chainid,resj_num,resj_name" << endl;

    fstream cont_out;
    MstUtils::openFile(cont_out, structure_name+"_drawcontacts.tsv",fstream::out);

    fstream ss_out;
    MstUtils::openFile(ss_out, structure_name+"_drawcontacts_ss.tsv",fstream::out);

    fstream sb_out;
    MstUtils::openFile(sb_out, structure_name+"_drawcontacts_sb.tsv",fstream::out);

    fstream bs_out;
    MstUtils::openFile(bs_out, structure_name+"_drawcontacts_bs.tsv",fstream::out);

    fstream bb_out;
    MstUtils::openFile(bb_out, structure_name+"_drawcontacts_bb.tsv",fstream::out);
    
    int count = 0;
    for (Residue* Ri : S.getResidues()) {
        set<Residue*> contactingResiduesSS = vdwC.getInteractingRes(Ri,vdwContacts::vdwContactType::SIDECHAIN);
        for (Residue* Rj : contactingResiduesSS) {
            out << S.getName() << ",";
            out << Ri->getChainID() << "," << Ri->getNum() << "," << Ri->getName() << ",";
            out << Rj->getChainID() << "," << Rj->getNum() << "," << Rj->getName() << ",";
            out << endl;

            ss_out << Ri->getChainID() << Ri->getNum() << "\t";
            ss_out << Rj->getChainID() << Rj->getNum() << "\t";
            ss_out << 1;
            ss_out << endl;

            count++;
        }

        set<Residue*> contactingResiduesSB = vdwC.getInteractingRes(Ri,vdwContacts::vdwContactType::SIDECHAIN_BACKBONE);
        for (Residue* Rj : contactingResiduesSB) {
            out << S.getName() << ",";
            out << Ri->getChainID() << "," << Ri->getNum() << "," << Ri->getName() << ",";
            out << Rj->getChainID() << "," << Rj->getNum() << "," << Rj->getName() << ",";
            out << endl;

            sb_out << Ri->getChainID() << Ri->getNum() << "\t";
            sb_out << Rj->getChainID() << Rj->getNum() << "\t";
            sb_out << 1;
            sb_out << endl;

            count++;
        }

        set<Residue*> contactingResiduesBS = vdwC.getInteractingRes(Ri,vdwContacts::vdwContactType::BACKBONE_SIDECHAIN);
        for (Residue* Rj : contactingResiduesBS) {
            out << S.getName() << ",";
            out << Ri->getChainID() << "," << Ri->getNum() << "," << Ri->getName() << ",";
            out << Rj->getChainID() << "," << Rj->getNum() << "," << Rj->getName() << ",";
            out << endl;

            bs_out << Ri->getChainID() << Ri->getNum() << "\t";
            bs_out << Rj->getChainID() << Rj->getNum() << "\t";
            bs_out << 1;
            bs_out << endl;

            count++;
        }
        set<Residue*> contactingResiduesBB = vdwC.getInteractingRes(Ri,vdwContacts::vdwContactType::BACKBONE);
        for (Residue* Rj : contactingResiduesBB) {
            out << S.getName() << ",";
            out << Ri->getChainID() << "," << Ri->getNum() << "," << Ri->getName() << ",";
            out << Rj->getChainID() << "," << Rj->getNum() << "," << Rj->getName() << ",";
            out << endl;

            bb_out << Ri->getChainID() << Ri->getNum() << "\t";
            bb_out << Rj->getChainID() << Rj->getNum() << "\t";
            bb_out << 1;
            bb_out << endl;

            count++;
        }
    }
    cout << "There are " << count << " vdw contacts" << endl;
    out.close();
    cont_out.close();

    return 0;
}
