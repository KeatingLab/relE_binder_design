#ifndef _RESIDUECONTACTS_H
#define _RESIDUECONTACTS_H

#include "msttypes.h"
#include <set>

// Craig's vdwRadii functions
enum atomInteractionType { ATOMCLASH = 0, ATOMCONTACT = 1, NOINTERACTION = 2};

class checkVDWRadii {
public:
    checkVDWRadii(bool verbose = true) {initConstants();};

    // /*
    // Renames each atom to a name in 'radii', or discards the atom (e.g. if hydrogen)
    // @param S The structure that is to processed
    // @param strict If true, will terminate if unknown atom encountered
    // */
    // void processStructure(Structure& S, bool strict = true);

    mstreal maxRadius();
    mstreal maxSumRadii();
    mstreal getRadius(const Atom& a1, bool strict = false);
    mstreal sumRadii(const Atom& a1, const Atom& a2);
    bool clash(const Atom& a1, const Atom& a2, mstreal lb = 0.7);
    bool contact(const Atom& a1, const Atom& a2, mstreal lb = 0.7, mstreal ub = 1.0);
    bool independent(const Atom& a1, const Atom& a2, mstreal ub = 1.0);
    atomInteractionType interactionType(const Atom& a1, const Atom& a2, mstreal lb = 0.7, mstreal ub = 1.0);
private:
    bool initConstants();
    string getName(const Atom* a, bool strict = false);
    
    map<string, mstreal> radii;
    set<string> ignore_atoms;
    mstreal max_radius;
    bool verbose;

};

class vdwContacts {
    /* A light-weight class that can identify pairs of residues that make VDW interactions (through
    their sidechains or otherwise).
     
     NOTE:
     Assumes all atom names in the structure are recognizable
     
     */
public:
    vdwContacts(Structure* _S);
    
    set<Residue*> getInteractingRes(Residue* R);

    map<int,set<int> > getAllInteractingRes();
    
protected:
    
private:
    Structure* S;
    vector<Atom*> allAtoms;
    ProximitySearch allAtomsPS;
    
    checkVDWRadii vdwR;
        
    map<Residue*,set<Residue*>> interacting;
    
    int ignoreDistance = 8; //interacting residues this close in the chain will be ignored
    mstreal adjustment = 1.0; //makes the interaction check more/less permissive
    bool strictAtomName = false;
};

#endif