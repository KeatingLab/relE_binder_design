#include "mstsystem.h"
#include "mstrotlib.h"
#include "mstoptions.h"

#include "residueframe.h"

int main(int argc, char *argv[]) {
    MstOptions op;
    op.setTitle("Reads protein structures and gets the Cb vector spherical coordinates");
    op.addOption("pdbList", "List of PDB structures",true);
    op.setOptions(argc,argv);

    string pdbListPath = op.getString("pdbList");
    vector<string> pdbList = MstUtils::fileToArray(pdbListPath);

    fstream out;
    MstUtils::openFile(out,"structuresCbSphericalCoordinates.csv",fstream::out);
    out << "structure,chainID,resNum,resName,radius,polarAngle,azimuthalAngle" << endl;

    TransformFactory tf;
    for (string structureName : pdbList) {
        Structure S(structureName);
        string name = MstSys::splitPath(S.getName(),1);
        for (Residue* R : S.getResidues()) {
            if (!RotamerLibrary::hasFullBackbone(R)) {
                cout << "Skipping residue: " << R->getChainID() << R->getNum() << " in structure " << R->getParent()->getStructure()->getName() << endl;
                continue;
            }
            Atom* Cb = R->findAtom("CB",false);
            Atom* Ca = R->findAtom("CA",false);
            if ((Cb == NULL) || (Ca == NULL)) {
                cout << "Skipping residue type: " << R->getName() << endl;
                continue;
            }
            CartesianPoint CaCb = Cb->getCoor() - Ca->getCoor();

            residueFrame rF(R);
            rF.setO(CartesianPoint(0,0,0));
            Transform localToGlobal = tf.switchFrames(Frame(),rF);
            localToGlobal.apply(CaCb);

            mstreal radius, polarAngle, azimuthalAngle = 0.0;
            CaCb.convertToSphericalCoordinates(radius,polarAngle,azimuthalAngle);

            out << name << "," << R->getChainID() << ",";
            out << R->getNum() << "," << R->getName() << ",";
            out << radius << "," << polarAngle << ",";
            out << azimuthalAngle << endl;
        }

    }
    return 0;
}