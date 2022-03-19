#include "mstoptions.h"
#include "mstrotlib.h"
#include "mstsystem.h"

#include "residuecontact.h"
#include "residueframe.h"

int main(int argc, char *argv[]) {
    MstOptions op;
    op.setTitle("Writes contact data from structures. Finds all potential contacts within certain distance and reports their cosine angle and whether they are real contacts");
    op.addOption("pdbList", "List of PDB structures",true);
    op.addOption("omitNonContacting", "If provided, will report real VDW contacts (so all potential contacts with no VDW contact will be omitted)");
    op.setOptions(argc,argv);

    bool omitNonContact = op.isGiven("omitNonContacting");
    string pdbListPath = op.getString("pdbList");
    vector<string> pdbList = MstUtils::fileToArray(pdbListPath);

    fstream out;
    string fileName = MstSys::splitPath(pdbListPath,1) + ".csv";
    MstUtils::openFile(out,fileName,fstream::out);

    out << "structure,resIchainID,resInum,resIname,resJchainID,resJnum,resJname,CaDistance,normalizedCbdistance,RiCbOrientation,RjCbOrientation,vdwContactSS,vdwContactSB,vdwContactBS,vdwContactBB" << endl;

    for (string pdbPath : pdbList) {
        Structure S(pdbPath);
        S.setName(MstSys::splitPath(S.getName(),1));

        // remove all residues with incomplete backbone
        vector<Residue*> validRes;
        for (Residue* R : S.getResidues()) {
            if (!RotamerLibrary::hasFullBackbone(R)) {
                cout << "Skipping residue: " << R->getChainID() << R->getNum() << " in structure " << R->getParent()->getStructure()->getName() << endl;
                continue;
            }
            validRes.push_back(R);
        }

        vdwContacts vdwC(validRes);
        map<int,set<int>> vdwSSConts = vdwC.getAllInteractingRes(vdwContacts::vdwContactType::SIDECHAIN);
        map<int,set<int>> vdwSBConts = vdwC.getAllInteractingRes(vdwContacts::vdwContactType::SIDECHAIN_BACKBONE);
        map<int,set<int>> vdwBSConts = vdwC.getAllInteractingRes(vdwContacts::vdwContactType::BACKBONE_SIDECHAIN);
        map<int,set<int>> vdwBBConts = vdwC.getAllInteractingRes(vdwContacts::vdwContactType::BACKBONE);
        cout << "Found " << vdwSSConts.size() + vdwSBConts.size() + vdwBSConts.size() + vdwBBConts.size() << " vdw contacts" << endl;

        potentialContacts pC(validRes,validRes);
        vector<pair<Residue*,Residue*>> pContacts = pC.getContacts(true);
        cout << "Found " << pContacts.size() << " distance based contacts" << endl;

        for (auto contact : pContacts) {
            Residue* targetRes = contact.first;
            Residue* binderRes = contact.second;

            bool isvdwSSCont = false, isvdwSBCont = false, isvdwBSCont = false, isvdwBBCont = false;
            map<int,set<int>>::iterator it = vdwSSConts.find(targetRes->getResidueIndex());
            if ((it != vdwSSConts.end())&&(it->second.find(binderRes->getResidueIndex()) != it->second.end())) {
                isvdwSSCont = true;
            }
            map<int,set<int>>::iterator it2 = vdwSBConts.find(targetRes->getResidueIndex());
            if ((it2 != vdwSBConts.end())&&(it2->second.find(binderRes->getResidueIndex()) != it2->second.end())) {
                isvdwSBCont = true;
            }
            map<int,set<int>>::iterator it3 = vdwBSConts.find(targetRes->getResidueIndex());
            if ((it3 != vdwBSConts.end())&&(it3->second.find(binderRes->getResidueIndex()) != it3->second.end())) {
                isvdwBSCont = true;
            }
            map<int,set<int>>::iterator it4 = vdwBBConts.find(targetRes->getResidueIndex());
            if ((it4 != vdwBBConts.end())&&(it4->second.find(binderRes->getResidueIndex()) != it4->second.end())) {
                isvdwBBCont = true;
            }

            if (omitNonContact & (!isvdwSSCont & !isvdwSBCont & !isvdwBSCont & !isvdwBBCont)) continue;

            mstreal CaDistance = targetRes->findAtom("CA")->getCoor().distance(binderRes->findAtom("CA"));
            mstreal normCbDistance = pC.getNormalizedCbDistance(targetRes,binderRes);
            mstreal RiCbOrientation = pC.getCaCbtoRiCaRjCaAngle(targetRes,binderRes);
            mstreal RjCbOrientation = pC.getCaCbtoRiCaRjCaAngle(binderRes,targetRes);

            out << S.getName() << "," << targetRes->getChainID() << "," << targetRes->getNum() << "," << targetRes->getName() << ",";
            out << binderRes->getChainID() << "," << binderRes->getNum() << "," << binderRes->getName() << ",";
            out << CaDistance << "," << normCbDistance << "," << RiCbOrientation << "," << RjCbOrientation << ",";
            out << isvdwSSCont << "," << isvdwSBCont << "," << isvdwBSCont << "," << isvdwBBCont << endl;
        }
    }
    out.close();

    return 0;
}