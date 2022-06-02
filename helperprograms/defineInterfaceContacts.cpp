#include "msttypes.h"
#include "mstoptions.h"

#include "alignframes.h"
#include "generateseeds.h"
#include "scoreinterface.h"
#include "residuecontact.h"
#include "residueframe.h"

#include <chrono>

int main(int argc, char *argv[]) {
    MstOptions op;
    op.setTitle("Defines contacts (either heavy-atom, or potential) between residues in the binder chain and residues in all other chains.");
    op.addOption("contactData","Path to a JSON file describing probability density of contacts",true);
    op.addOption("complexPDB","A PDB file defining the complex between the protein target and designed binder. The target must have a defined sequence",false);
    op.addOption("binderChainID","The ID of the binder chain",true);
    op.addOption("vdwContacts","If provided, will define the interface using VDW contacts, otherwise will use CB definition. Not compatible with --targetPDB mode",false);
    op.setOptions(argc,argv);

    string contactData = op.getString("contactData");
    string complexPDBPath = op.getString("complexPDB","");
    string binderChainID = op.getString("binderChainID");
    bool defineVDWContacts = op.isGiven("vdwContacts");

    string fileType = (defineVDWContacts) ? "vdw" : "potential";

    fstream cont_out;
    MstUtils::openFile(cont_out,fileType+"_AllContacts.tsv",fstream::out);
    // cont_out << "binderChainID,binderResNum,targetChainID,targetResNum" << endl;
    // No header, file intended to be read by contact drawing pymol script

    fstream clash_out;
    MstUtils::openFile(clash_out,fileType+"_AllClashes.tsv",fstream::out);

    fstream cont_byres_out;
    MstUtils::openFile(cont_byres_out,fileType+"_BinderResContacts.csv",fstream::out);
    cont_byres_out << "binderChainID,binderResNum,numContacts" << endl;

    Structure complexPDB(complexPDBPath);
    vector<Chain*> binderChains;
    vector<Chain*> targetChains;
    for (int i = 0; i < complexPDB.chainSize(); i++) {
        Chain* C = &complexPDB.getChain(i);
        if (C->getID() == binderChainID) binderChains.push_back(C);
        else targetChains.push_back(C);
    }

    if (defineVDWContacts) {
        vdwContacts conts(binderChains,targetChains);
        vector<pair<Residue*,Residue*>> contPairs = conts.getInteractingResPairs();
        for (auto cont : contPairs) {
            cont_out << cont.first->getChainID() << cont.first->getNum() << "\t";
            cont_out << cont.second->getChainID() << cont.second->getNum() << "\t";
            cont_out << 1.0 << endl;
        }

        vector<Residue*> binderResidues = complexPDB.getChainByID(binderChainID)->getResidues();
        for (Residue* binderRes : binderResidues) {
            set<Residue*> contactingRes = conts.getInteractingRes(binderRes);
            cont_byres_out << binderRes->getChainID() << "," << binderRes->getNum() << ",";
            cont_byres_out << contactingRes.size() << endl;
        }

    } else {
        vector<Residue*> binderRes, targetRes;
        binderRes = complexPDB.getChainByID(binderChainID)->getResidues();
        for (Chain* C : targetChains) {
            vector<Residue*> chainRes = C->getResidues();
            targetRes.insert(targetRes.end(),chainRes.begin(),chainRes.end());
        }
        potentialContacts potConts(targetRes,binderRes);
        potConts.load2DProbabilityDensities(contactData);

        vector<pair<Residue*,Residue*>> contPairs = potConts.getContacts();
        for (auto cont : contPairs) {
            cont_out << cont.first->getChainID() << cont.first->getNum() << "\t";
            cont_out << cont.second->getChainID() << cont.second->getNum() << "\t";
            cont_out << 1.0 << endl;
        }

        vector<pair<Residue*,Residue*>> clashPairs = potConts.getNonDesignableContacts();
        cout << "Found " << clashPairs.size() << " pairs of non-designable residues" << endl;
        for (auto cont : clashPairs) {
            clash_out << cont.first->getChainID() << cont.first->getNum() << "\t";
            clash_out << cont.second->getChainID() << cont.second->getNum() << "\t";
            clash_out << 1.0 << endl;
        }

        vector<Residue*> binderResidues = complexPDB.getChainByID(binderChainID)->getResidues();
        for (Residue* binderRes : binderResidues) {
            set<Residue*> contactingRes = potConts.getContactsWithResidue(binderRes);
            cont_byres_out << binderRes->getChainID() << "," << binderRes->getNum() << ",";
            cont_byres_out << contactingRes.size() << endl;
        }
    }
    cont_out.close();
    cont_byres_out.close();

    cout << "Done" << endl;

    return 0;
}