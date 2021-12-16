#include "residueframe.h"

residueFrame::residueFrame(Residue* R) : parent(R) {defineFrame(R);}

void residueFrame::defineFrame(Residue* R) {
    /*
    The position is defined as the 3D coordinate of the alpha-carbon, relative to the origin
    */
    bool strict = true;
    position = R->findAtom("CA",strict);

    /*
    Following Ingraham et. al 2019, the orientation is defined as three basis vectors: [b, n, b x n]
    b = the negative bisector of the angle between N - Ca and C - Ca
    n = the unit vector normal to the plane formed between N - Ca and C - Ca
    */

    CartesianPoint N = R->findAtom("N",strict);
    CartesianPoint Ca = position;
    CartesianPoint C = R->findAtom("C",strict);

    CartesianPoint N_to_Ca = (Ca-N)/(Ca-N).norm();
    CartesianPoint Ca_to_C = (C-Ca)/(C-Ca).norm();

    CartesianPoint n = N_to_Ca.cross(Ca_to_C)/(N_to_Ca.cross(Ca_to_C)).norm(); //n is x
    CartesianPoint b = (N_to_Ca-Ca_to_C)/(N_to_Ca-Ca_to_C).norm(); //b is y
    CartesianPoint u = n.cross(b)/(n.cross(b).norm());  //u is z

    mstreal x_component = u.getX() + b.getX() + n.getX();
    mstreal y_component = u.getY() + b.getY() + n.getY();
    mstreal z_component = u.getZ() + b.getZ() + n.getZ();
    bool unit_sphere = true;
    orientation.constructSphericalCoordinatesFromXYZ(x_component,y_component,z_component,unit_sphere);

    constructFrame(position,n,b,u);
}

residueFrame residueFrame::frameRelativeToOther(residueFrame* other) {
    // find position of other residue, relative to self
    CartesianPoint relPos = other->getPosition() - position;

    // find orientation of other residue, relative to self
    sphericalCoordinate relOr = other->getOrientation() - orientation;

    return residueFrame(relPos,relOr);
}

void augmentedStructure::writeToFile(string path_prefix) {
    // write each residue frame out such that the basis vector can be interpreted by cgo_basis.py script
    fstream out;
    MstUtils::openFile(out, path_prefix+"_frames.tsv",fstream::out);
    
    for (residueFrame& rF : frames) {
        Residue* R = rF.getParent();
        out << R->getChainID() << R->getNum() << "\t";
        out << rF.getPosition() << "\t";
        out << rF.getUPos() << "\t";
        out << rF.getBPos() << "\t";
        out << rF.getNPos();
        out << endl;
    }
    out.close();
}