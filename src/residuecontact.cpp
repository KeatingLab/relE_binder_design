#include "residuecontact.h"

// void checkVDWRadii::processStructure(Structure& S, bool strict) {
//     set<Atom*> toDelete;
//     for (Residue* R : S.getResidues()) {
//         vector<Atom*> atoms = R->getAtoms();
//         for (Atom* A : atoms) {
//             string name = getName(A);
//             if (name.empty()) toDelete.insert(A);
//             A->setName(name);
//         }
//     }
// }

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

set<Residue*> vdwContacts::getInteractingRes(Residue* R) {
    // Check if map is already populated
    if (interacting.count(R) == 1) return interacting[R];
    // If not, compute the contacts anew
    vector<Atom*> resAtoms = R->getAtoms();
    set<Residue*> contacts;
    for (Atom* A : resAtoms) {
        mstreal A_maxRadius_sum = vdwR.getRadius(A,strictAtomName) + vdwR.maxRadius();
        vector<int> atomsToCheck = allAtomsPS.getPointsWithin(A->getCoor(), 0, A_maxRadius_sum);
        for (int atomIdx : atomsToCheck) {
            Atom* A_other = allAtoms[atomIdx];
            Residue* R_other = A_other->getResidue();
            
            // Check if should ignore, since in the same residue
            if (R == R_other) {
                continue;
            }
            
            // Check if should ignore, since too close in chain
            if (R->getParent()==R_other->getParent()) {
                int R_i = R->getResidueIndexInChain();
                int R_j = R_other->getResidueIndexInChain();
                int distance = abs(R_i-R_j);
                if (distance <= ignoreDistance) {
                    continue;
                }
            }
            
            // Check if contacting, given the Aj radius
            atomInteractionType intType = vdwR.interactionType(A,A_other);
            if (intType == ATOMCONTACT) {
                contacts.insert(R_other);
            }
        }
    }
    return contacts;
}

map<int,set<int> > vdwContacts::getAllInteractingRes(bool verbose) {
    map<int,set<int> > allInteracting;
    int numContacts = 0;
    for (Residue* Ri : S.getResidues()) {
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