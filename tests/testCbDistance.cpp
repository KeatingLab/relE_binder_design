#include "mstoptions.h"
#include "mstrotlib.h"
#include "mstsystem.h"

#include "residuecontact.h"
#include "residueframe.h"

int main(int argc, char *argv[]) {
    MstOptions op;
    op.setTitle("Computes normalized Cb distance between residues and writes a PDB file with inferred Cb coordinates");
    op.addOption("pdb", "structure to compute contacts for",true);
    op.setOptions(argc,argv);

    string pdbPath = op.getString("pdb");

    Structure inputStructure(pdbPath);
    Structure bbStructure(RotamerLibrary::getBackbone(inputStructure));
    cout << "extracted backbone with " << bbStructure.residueSize() << " residues" << endl;

    potentialContacts pC(bbStructure.getResidues(),bbStructure.getResidues(),true,false);
    vector<pair<Residue*,Residue*>> pContacts = pC.getContacts(true);
    cout << "Found " << pContacts.size() << " distance based contacts" << endl;

    for (auto contact : pContacts) {
        Residue* targetRes = contact.first;
        Residue* binderRes = contact.second;

        // bool isVDWcont = false;
        // map<int,set<int>>::iterator it = vdwConts.find(targetRes->getResidueIndex());
        // if ((it != vdwConts.end())&&(it->second.find(binderRes->getResidueIndex()) != it->second.end())) {
        //     isVDWcont = true;
        // }

        mstreal CaDistance = targetRes->findAtom("CA")->getCoor().distance(binderRes->findAtom("CA"));
        mstreal normCbDistance = pC.getNormalizedCbDistance(targetRes,binderRes);

        cout << targetRes->getChainID() << targetRes->getNum() << "\t";
        cout << binderRes->getChainID() << binderRes->getNum() << "\t";
        cout << CaDistance << "\t" << normCbDistance << endl;
    }

    // add Cb to bb structure
    for (Residue* R: bbStructure.getResidues()) {
        CartesianPoint CbPoint = pC.getCbFromRes(R);

        Atom* Cb = new Atom();
        Cb->setX(CbPoint.getX());
        Cb->setY(CbPoint.getY());
        Cb->setZ(CbPoint.getZ());
        Cb->setName("CB");

        R->appendAtom(Cb);
    }
    bbStructure.writePDB("bbStructureWCb.pdb");

    return 0;
}