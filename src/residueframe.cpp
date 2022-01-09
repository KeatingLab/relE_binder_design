#include "residueframe.h"

/* --- --- --- --- --- residueFrame --- --- --- --- --- */

residueFrame* residueFrame::frameRelativeToOther(const residueFrame& other) {
    // find the transformation between other and the global reference frame
    Transform other_to_ref = TransformFactory::switchFrames(Frame(), other);

    // apply the transformation to this frame and return
    residueFrame* result = new residueFrame(*this);
    other_to_ref.apply(*result);
    return result;
}

void residueFrame::writeToFile(string name, fstream& out) {
    out << name << "\t";
    out << this->getO() << "\t";
    out << this->getXPos() << "\t";
    out << this->getYPos() << "\t";
    out << this->getZPos();
    out << endl;
}

void residueFrame::defineFrame(Residue* R) {
    /*
    Following Ingraham et. al 2019, the orientation is defined as three basis vectors: [b, n, b x n]
    b = the negative bisector of the angle between N - Ca and C - Ca
    n = the unit vector normal to the plane formed between N - Ca and C - Ca
    */

    bool strict = true;
    CartesianPoint N = R->findAtom("N",strict);
    CartesianPoint Ca = R->findAtom("CA",strict);
    CartesianPoint C = R->findAtom("C",strict);

    CartesianPoint N_to_Ca = (Ca-N);
    CartesianPoint Ca_to_C = (C-Ca);

    CartesianPoint z = Ca_to_C.cross(N_to_Ca); //normal to residue plane (n)
    CartesianPoint y = (N_to_Ca-Ca_to_C); // (b)
    CartesianPoint x = y.cross(z);  //(u)

    constructFrame(Ca,x,y,z);
}

/* --- --- --- --- --- augmentedStructure --- --- --- --- --- */

void augmentedStructure::writeToFile(string path_prefix) {
    // write each residue frame out such that the basis vector can be interpreted by cgo_basis.py script
    fstream out;
    MstUtils::openFile(out, path_prefix+"_frames.tsv",fstream::out);
    
    for (residueFrame& rF : frames) {
        Residue* R = rF.getParent();
        string name = R->getChainID() + MstUtils::toString(R->getNum());
        rF.writeToFile(name,out);
    }
    out.close();
}