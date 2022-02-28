#include "mstoptions.h"
#include "mstrotlib.h"
#include "mstsystem.h"

#include "residuecontact.h"
#include "residueframe.h"

int main(int argc, char *argv[]) {
    MstOptions op;
    op.setTitle("Writes contact data from structures. Finds all potential contacts within certain distance and reports their cosine angle and whether they are real contacts");
    op.addOption("pdbList", "List of PDB structures",true);
    op.addOption("writeCb","write structure with Cb atoms");
    op.setOptions(argc,argv);

    string pdbListPath = op.getString("pdbList");
    vector<string> pdbList = MstUtils::fileToArray(pdbListPath);

    fstream out;
    string fileName = MstSys::splitPath(pdbListPath,1) + ".csv";
    MstUtils::openFile(out,fileName,fstream::out);

    out << "structure,resIchainID,resInum,resIname,resJchainID,resJnum,resJname,CaDistance,normalizedCbdistance,vdwContact" << endl;

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
        map<int,set<int>> vdwConts = vdwC.getAllInteractingRes();
        cout << "Found " << vdwConts.size() << " vdw contacts" << endl;

        potentialContacts pC(validRes,validRes,true,false);
        vector<pair<Residue*,Residue*>> pContacts = pC.getContacts(true);
        cout << "Found " << pContacts.size() << " distance based contacts" << endl;

        for (auto contact : pContacts) {
            Residue* targetRes = contact.first;
            Residue* binderRes = contact.second;

            bool isVDWcont = false;
            map<int,set<int>>::iterator it = vdwConts.find(targetRes->getResidueIndex());
            if ((it != vdwConts.end())&&(it->second.find(binderRes->getResidueIndex()) != it->second.end())) {
                isVDWcont = true;
            }

            mstreal CaDistance = targetRes->findAtom("CA")->getCoor().distance(binderRes->findAtom("CA"));
            mstreal normCbDistance = pC.getNormalizedCbDistance(targetRes,binderRes);

            out << S.getName() << "," << targetRes->getChainID() << "," << targetRes->getNum() << "," << targetRes->getName() << ",";
            out << binderRes->getChainID() << "," << binderRes->getNum() << "," << binderRes->getName() << ",";
            out << CaDistance << "," << normCbDistance << "," << isVDWcont << endl;
        }
    }
    out.close();

    return 0;
}