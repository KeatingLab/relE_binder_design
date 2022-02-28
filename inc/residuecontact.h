#ifndef _RESIDUECONTACTS_H
#define _RESIDUECONTACTS_H

#include <set>

#include "mstsequence.h"
#include "msttypes.h"

#include "residueframe.h"
#include "freesasaext.h"
#include "utilities.h"

// Craig's vdwRadii functions
enum atomInteractionType { ATOMCLASH = 0, ATOMCONTACT = 1, NOINTERACTION = 2};

class checkVDWRadii {
public:
    checkVDWRadii(bool verbose = true) {initConstants();};

    mstreal maxRadius();
    mstreal maxSumRadii();
    mstreal getRadius(const Atom& a1, bool strict = false);
    mstreal getRadius(const string& resName, const string& atomName, bool strict = false);
    mstreal sumRadii(const Atom& a1, const Atom& a2);
    bool clash(const Atom& a1, const Atom& a2, mstreal lb = 0.7);
    bool contact(const Atom& a1, const Atom& a2, mstreal lb = 0.7, mstreal ub = 1.0);
    bool independent(const Atom& a1, const Atom& a2, mstreal ub = 1.0);
    atomInteractionType interactionType(const Atom& a1, const Atom& a2, mstreal lb = 0.7, mstreal ub = 1.0);
private:
    bool initConstants();
    string getName(string name, bool strict = false);
    
    // map<string, mstreal> radii;
    map<string,map<string,mstreal>> radii;
    set<string> ignore_atoms;
    mstreal max_radius = 0.0;
    bool verbose;

};

class vdwContacts {
    /* A light-weight class that can identify pairs of residues that make VDW interactions (through
    their sidechains or other atoms).
    
    This class has two uses: 
    1) The user wants to identify all contacts within a structure
    2) The user wants to identify all contacts between two (non-overlapping) sets of residues

    The use case is determined by the constructor that is selected. Admittedly it is a little weird to handle the two
    situations this way, but it will work for now.
     
     */
public:
    vdwContacts(vector<Residue*> S_res);
    vdwContacts(vector<Chain*> resIChains, vector<Chain*> resJChains);
    vdwContacts(vector<Residue*> resIVec, vector<Residue*> resJVec);
    
    set<Residue*> getInteractingRes(Residue* Ri);

    vector<pair<Residue*,Residue*>> getInteractingRes();

    map<int,set<int> > getAllInteractingRes(bool verbose = false);
    
protected:
    void setResidues(vector<Residue*> resIVec, vector<Residue*> resJVec);
    void preparePS(vector<Residue*> toAdd = {});

    bool isContact(Residue* Ri, Residue* Rj);
    
private:
    vector<Atom*> allAtoms;
    ProximitySearch allAtomsPS;
    
    checkVDWRadii vdwR;
        
    set<Residue*> resISet;
    set<Residue*> resJSet;
    map<Residue*,set<Residue*>> interacting;
    
    int ignoreDistance = 8; //interacting residues this close in the chain will be ignored
    mstreal adjustment = 1.0; //makes the interaction check more/less permissive
    bool strictAtomName = false;
};

struct potentialContactsParams {
    vector<mstreal> CAdistanceToCheck; // each residue type has a different distance
    vector<mstreal> decayRate;
};

class potentialContacts {
    /* Identifies pairs of residues that are positioned such that they could form a contact
     
     This class is designed exclusively for the scenario where we are considering a target protein and binder.
     In this case, the residue identity of the target protein residues is fixed.
     */

    public:
        potentialContacts(vector<Residue*> _targetResidues, vector<Residue*> _binderResidues, bool sameChain = false, bool strict = true);

        vector<pair<Residue*,Residue*>> getContacts(bool simpleDistanceCheck = false);

        CartesianPoint getCbFromRes(Residue* R);
        mstreal getNormalizedCbDistance(Residue* Ri, Residue* Rj);

    protected:
        bool ignoreContact(Residue* Ri, Residue* Rj);

    private:
        vector<Residue*> targetResidues;
        vector<Residue*> binderResidues;

        vector<CartesianPoint> targetResCb;
        vector<CartesianPoint> binderResCb;

        vector<Atom*> targetCA;
        vector<Atom*> binderCA;

        // estimated from 1000 randomly sampled structures (using tests/computeCbParams.cpp)
        mstreal radius = 1.5;
        mstreal polarAngle = 37.794;
        mstreal azimuthalAngle = 87.948;

        map<Residue*,res_t> targetResidueAAIdentity;
        mstreal defaultDistanceToCheck = 10.0;

        bool sameChain; // if true, will check whether nearby residues are in the same chain
        int ignoreDistance = 8; //interacting residues this close in the chain will be ignored
};

#endif