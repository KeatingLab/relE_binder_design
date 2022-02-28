#include "utilities.h"

/* --- --- --- --- --- SeqToolsExtension --- --- --- --- --- */

// statics must be defined in the cpp file
set<string> SeqToolsExtension::aaNames;
map<string,string> SeqToolsExtension::equivalentAA;
bool SeqToolsExtensionInitialized = SeqToolsExtension::initConstants();

bool SeqToolsExtension::initConstants() {
    // for now, will aggregate data into the 20 canonical amino acids
    aaNames = {"ALA","ARG","ASN","ASP",
                "CYS","GLN","GLU","GLY",
                "HIS","ILE","LEU","LYS",
                "MET","PHE","PRO","SER",
                "THR","TRP","TYR","VAL"};

    // in some cases, we may want to consider these (identical/similar) residues
    //// equivalent
    equivalentAA["HSD"] = "HIS";
    equivalentAA["HSE"] = "HIS";
    equivalentAA["HSC"] = "HIS";
    equivalentAA["HSP"] = "HIS";

    //// similar
    equivalentAA["MSE"] = "MET";
    equivalentAA["SEC"] = "CYS";
    return true;
}

set<string> SeqToolsExtension::getAANames() {
    return aaNames;
}

string SeqToolsExtension::findEquivalentResidueInAlphabet(string aa3, bool strict) {
    bool found = (equivalentAA.find(aa3) != equivalentAA.end());
    if (!found) {
        if (strict) MstUtils::error("No residue equivalent to "+aa3,"SeqToolsExtension::findEquivalentResidueInAlphabet");
        else return "";
    }
    return equivalentAA[aa3];
}

bool SeqToolsExtension::AAinSet(string aa3) {
    return (aaNames.find(aa3) != aaNames.end());
}

int SeqToolsExtension::numAA() {
    return aaNames.size();
}

/* --- --- --- --- --- MiscTools --- --- --- --- --- */

MiscTools::alignment MiscTools::bestRMSD(Chain* Ci, Chain* Cj, int kmerL) {
        if (kmerL > min(Ci->residueSize(),Cj->residueSize())) kmerL = min(Ci->residueSize(),Cj->residueSize());
    alignment bestAlignment;
    mstreal newRMSD;
    vector<Atom*> CiAtoms = getBackboneAtoms(Ci);
    vector<Atom*> CjAtoms = getBackboneAtoms(Cj);
    for (int i = 0; i < Ci->residueSize() - kmerL + 1; i++) {
        vector<Atom*> CiAtomsWindowi(CiAtoms.begin()+(4*i),CiAtoms.begin()+(4*(i+kmerL)));
        for (int j = 0; j < Cj->residueSize() - kmerL + 1; j++) {
            vector<Atom*> CjAtomsWindowj(CjAtoms.begin()+(4*j),CjAtoms.begin()+(4*(j+kmerL)));
            newRMSD = RMSDCalculator::rmsd(CiAtomsWindowi,CjAtomsWindowj);
            if (newRMSD < bestAlignment.rmsd) {
                bestAlignment.rmsd = newRMSD;
                bestAlignment.CiResIdx = i;
                bestAlignment.CjResIdx = j;
                bestAlignment.length = kmerL;
            }
        }
    }
    return bestAlignment;
}

void MiscTools::extractBackboneFromStructure(const Structure& source, Structure& destination) {
    if (destination.atomSize() > 0) {
        cout << "Clearing all atoms from " << destination.getName() << endl;
        destination.reset();
    }
    for (Residue* R : source.getResidues()) {
        if (!RotamerLibrary::hasFullBackbone(R)) MstUtils::error("Source structure missing backbone atoms.","MiscTools::extractBackboneFromStructure");
        vector<Atom*> RbbAtoms = RotamerLibrary::getBackbone(R);
        for (Atom* bbAtom : RbbAtoms) destination.addAtom(bbAtom);
    }
    destination.setName(source.getName());
}