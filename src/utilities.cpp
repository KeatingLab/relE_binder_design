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
    equivalentAA["TYS"] = "TYR";
    equivalentAA["PTR"] = "TYR";

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

/* --- --- --- --- --- twoDimHistogram --- --- --- --- --- */

void twoDimHistogram::setCounts(vector<vector<int>> _counts) {
    if (_counts.size() != dim1NBins) MstUtils::error("Wrong number of bins along dimension 1","twoDimHistogram::setCounts");
    counts.clear();
    normalizingConstant = 0;
    for (vector<int> row : _counts) {
        if (row.size() != dim2NBins)  MstUtils::error("Wrong number of bins along dimension 2","twoDimHistogram::setCounts");
        for (int val : row) normalizingConstant+=val;
        counts.push_back(row);
    }
}

int twoDimHistogram::getCounts(mstreal val1, mstreal val2) {
    if (counts.empty()) MstUtils::error("Counts not set!","twoDimHistogram::getCounts");
    int idx1, idx2;
    getIdx(val1,val2,idx1,idx2);
    return counts[idx1][idx2];
}

mstreal twoDimHistogram::getDensity(mstreal val1, mstreal val2) {
    if (counts.empty()) MstUtils::error("Counts not set!","twoDimHistogram::getCounts");
    int idx1, idx2;
    getIdx(val1,val2,idx1,idx2);
    return mstreal(counts[idx1][idx2]) / mstreal(normalizingConstant) / (dim1BinSize*dim2BinSize);
}

void twoDimHistogram::getIdx(mstreal val1, mstreal val2, int& idx1, int& idx2) {
    if (counts.empty()) MstUtils::error("Must set counts before attempting to retrieve values","twoDimHistogram::getIdx");
    if (val1 < 0) idx1 = 0;
    else idx1 = min(dim1NBins-1,int(floor((val1-dim1Min)/dim1BinSize)));
    if (val2 < 0) idx2 = 0;
    else idx2 = min(dim2NBins-1,int(floor((val2-dim2Min)/dim2BinSize)));
}

void twoDimHistogram::reportHistogram() {
    cout << "dim1NBins: " << dim1NBins << endl;
    cout << "dim2NBins: " << dim2NBins << endl;
    cout << "dim1Min: " << dim1Min << endl;
    cout << "dim1Max: " << dim1Max << endl;
    cout << "dim2Min: " << dim2Min << endl;
    cout << "dim2Max: " << dim2Max << endl;
    cout << "dim1BinSize: " << dim1BinSize << endl;
    cout << "dim2BinSize: " << dim2BinSize << endl;
    cout << "normalizingConstant: " << normalizingConstant << endl;
    for (int i = 0; i < counts.size(); i++) {
        for (int j = 0; j < counts.size(); j++) {
            cout << counts[i][j] << " ";
        }
        cout << endl;
    }
}