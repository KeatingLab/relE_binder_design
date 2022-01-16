#include "residuecontact.h"

/* --- --- --- --- --- checkVDWRadii --- --- --- --- --- */

mstreal checkVDWRadii::maxRadius() {
    return max_radius;
}

mstreal checkVDWRadii::maxSumRadii() {
    return checkVDWRadii::maxRadius() * 2;
}

mstreal checkVDWRadii::getRadius(const Atom& a, bool strict) {
    string atomNameInMap = getName(&a, strict);
    if (atomNameInMap.empty()) return 0.0;
    return radii[atomNameInMap];
}

mstreal checkVDWRadii::sumRadii(const Atom& a1, const Atom& a2) {
    return (getRadius(a1)+getRadius(a2));
}

bool checkVDWRadii::clash(const Atom& a1, const Atom& a2, double lb) {
    return (a1.distance(a2) < checkVDWRadii::sumRadii(a1, a2) * lb);
}

bool checkVDWRadii::contact(const Atom& a1, const Atom& a2, double lb, double ub) {
    double dist = a1.distance(a2);
    double s = checkVDWRadii::sumRadii(a1, a2);
    return (dist >= s * lb && dist < s * ub);
}

bool checkVDWRadii::independent(const Atom& a1, const Atom& a2, double ub) {
    return (a1.distance(a2) >= checkVDWRadii::sumRadii(a1, a2) * ub);
}

atomInteractionType checkVDWRadii::interactionType(const Atom& a1, const Atom& a2, double lb, double ub) {
    double dist = a1.distance(a2);
    double s = checkVDWRadii::sumRadii(a1, a2);
    if (dist < s * lb) {
        return ATOMCLASH;
    } else if (dist < s * ub) {
        return ATOMCONTACT;
    } else {
        return NOINTERACTION;
    }
}

bool checkVDWRadii::initConstants() {
    // https://pubs.rsc.org/en/content/articlelanding/2013/DT/c3dt50599e
    radii = {{"N",  1.60}, {"NT", 1.60}, {"CA", 2.365},
             {"C", 2.10},  {"CT", 2.10}, {"O", 1.60}, 
             {"CB", 2.2350}, {"S", 1.80}};
    ignore_atoms = {"H"};
    max_radius = 0;
    for (auto it : radii) if (it.second > max_radius) max_radius = it.second;
    // cout << "max radius: " << max_radius << endl;
    return true;
}

string checkVDWRadii::getName(const Atom* A, bool strict) {
    string atomName = A->getName();
    while (atomName.length() >= 1) {
        if (radii.count(atomName) > 0) {
            return atomName;
        } else {
            if (strict) break;
            atomName = atomName.substr(0, atomName.size() - 1);
        }
    }
    if (strict) MstUtils::error("VDW radius not defined for atom name"+atomName,"checkVDWRadii::getName");
    return "";
}

/* --- --- --- --- --- vdwContacts --- --- --- --- --- */

void vdwContacts::setResidues(vector<Chain*> resIChains, vector<Chain*> resJChains) {
    vector<Residue*> resIVec, resJVec;
    for (Chain* C : resIChains) {
        vector<Residue*> cRes = C->getResidues();
        resIVec.insert(resIVec.end(),cRes.begin(),cRes.end());
    }
    for (Chain* C : resJChains) {
        vector<Residue*> cRes = C->getResidues();
        resJVec.insert(resJVec.end(),cRes.begin(),cRes.end());
    }
    setResidues(resIVec,resJVec);
}

void vdwContacts::setResidues(vector<Residue*> resIVec, vector<Residue*> resJVec) {
    resISet = set<Residue*>(resIVec.begin(),resIVec.end());
    resJSet = set<Residue*>(resJVec.begin(),resJVec.end());
    cout << "set I has " << resISet.size() << " and set J has " << resJSet.size() << endl;
}

vdwContacts::vdwContacts(const Structure& _S) : S(_S) {
    // Construct proximity search object
    allAtoms = S.getAtoms();
    allAtomsPS = ProximitySearch(allAtoms,1.5);
    
    // Search each residue for other residues making VDW contacts
    for (Residue* R : S.getResidues()) {
        set<Residue*> contactedResidues = getInteractingRes(R);
        interacting[R] = contactedResidues;
    }
}

set<Residue*> vdwContacts::getInteractingRes(Residue* Ri) {
    vector<Atom*> resAtoms = Ri->getAtoms();
    set<Residue*> contacts;
    for (Atom* Ai : resAtoms) {
        mstreal A_maxRadius_sum = vdwR.getRadius(Ai,strictAtomName) + vdwR.maxRadius();
        vector<int> atomsToCheck = allAtomsPS.getPointsWithin(Ai->getCoor(), 0, A_maxRadius_sum);
        for (int atomIdx : atomsToCheck) {
            Atom* Aj = allAtoms[atomIdx];
            Residue* Rj = Aj->getResidue();
            
            if (!isContact(Ri,Rj)) continue;
            
            // Check if contacting, given the Aj radius
            atomInteractionType intType = vdwR.interactionType(Ai,Aj);
            if (intType == ATOMCONTACT) {
                contacts.insert(Rj);
            }
        }
    }
    return contacts;
}

vector<pair<Residue*,Residue*>> vdwContacts::getInteractingRes() {
    vector<pair<Residue*,Residue*>> result;
    int numContacts = 0;
    for (Residue* Ri : resISet) {
        set<Residue*> interactingRes;
        for (Residue* Rj : getInteractingRes(Ri)) {
            interactingRes.insert(Rj);
            numContacts++;
        }
        for (Residue* Rj : interactingRes) {
            result.emplace_back(pair<Residue*,Residue*>(Ri,Rj));
        }
    }
    // if (verbose) cout << "In total, selected residues have " << numContacts << " VDW contacts" << endl;
    return result;
}

map<int,set<int> > vdwContacts::getAllInteractingRes(bool verbose) {
    map<int,set<int> > allInteracting;
    int numContacts = 0;
    for (Residue* Ri : S.getResidues()) {
        // Check if map is already populated
        if (interacting.count(Ri) == 1) continue;
        // If not, compute the contacts anew
        set<int> interactingResIdx;
        for (Residue* Rj : getInteractingRes(Ri)) {
            interactingResIdx.insert(Rj->getResidueIndex());
            numContacts++;
        }
        allInteracting[Ri->getResidueIndex()] = interactingResIdx;
    }
    if (verbose) cout << "In total, structure has " << numContacts << " VDW contacts" << endl;
    return allInteracting;
}

bool vdwContacts::isContact(Residue* Ri, Residue* Rj) {
    // Check if should ignore, since in the same residue
    if (Ri == Rj) return false;

    // Check that both residues are in the sets under consideration
    if (resISet.find(Ri) == resISet.end()) return false;
    if (resJSet.find(Rj) == resJSet.end()) return false;

    // Check that residues are not too close in the chain
    if (Ri->getParent() == Ri->getParent()) {
        int Ri_pos = Ri->getResidueIndexInChain();
        int Rj_pos = Rj->getResidueIndexInChain();
        bool tooClose = (abs(Ri_pos-Rj_pos) <= ignoreDistance);
        if (tooClose) return false;
    }
    return true;
}