#ifndef _RESIDUECONTACTS_H
#define _RESIDUECONTACTS_H

#include <set>

#include "mstsequence.h"
#include "msttypes.h"

#include "residueframe.h"
#include "freesasaext.h"
#include "utilities.h"

#include "nlohmann/json.hpp"
using json = nlohmann::json;

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
    their sidechains or other atoms). This class has two uses: 
    1) The user wants to identify all contacts within a structure
    2) The user wants to identify all contacts between two (non-overlapping) sets of residues
    The use case is determined by the constructor that is selected. Admittedly it is a little weird to handle the two
    situations this way, but it will work for now.

    Note that the contact definitions are "low-resolution". There is no attempt to quantify whether an interaction is energetically
    favorable or anything like that, all a "contact" means is that given the atomic coordinates provided, two residues are
    in close enough proximity for their heavy atoms to directly influence one another. This is probably reasonable when looking
    at experimentally determined structures, but whenever a designed structure is being assessed, this defintion should be
    paired with some kind of potential that more carefully assesses whether the interaction has a favorable enthalpy.
     */
public:
    enum vdwContactType {NONE, ALL, SIDECHAIN, SIDECHAIN_BACKBONE, BACKBONE_SIDECHAIN, BACKBONE};

    vdwContacts(vector<Residue*> S_res);
    vdwContacts(vector<Chain*> resIChains, vector<Chain*> resJChains);
    vdwContacts(vector<Residue*> resIVec, vector<Residue*> resJVec);
    
    set<Residue*> getInteractingRes(Residue* Ri, vdwContactType contType = ALL);
    vector<pair<Residue*,Residue*>> getInteractingResPairs(vdwContactType contType = ALL);
    map<int,set<int>> getAllInteractingRes(vdwContactType contType = ALL, bool verbose = false);
    
protected:
    void setResidues(vector<Residue*> resIVec, vector<Residue*> resJVec);
    void preparePS(vector<Residue*> toAdd = {});
    void findAllContacts();

    bool isValid(Residue* Ri, Residue* Rj);
    vdwContactType classifyContact(Atom* Ai, Atom* Aj);
    
private:
    vector<Atom*> allAtoms;
    ProximitySearch allAtomsPS;
    
    checkVDWRadii vdwR;
        
    set<Residue*> resISet;
    set<Residue*> resJSet;
    map<Residue*,set<Residue*>> allInteracting; //union of interacting residue pairs
    map<Residue*,set<Residue*>> sidechainInteracting;
    map<Residue*,set<Residue*>> sidechainBackboneInteracting;
    map<Residue*,set<Residue*>> backboneSidechainInteracting;
    map<Residue*,set<Residue*>> backboneInteracting;
    
    int ignoreDistance = 8; //interacting residues this close in the chain will be ignored
    mstreal lowerBound = 0.7; //adjust the interaction check 
    mstreal upperBound = 1.2; //adjust the interaction check
    bool strictAtomName = false;

    set<string> bbAtomNames = {"N","CA","C","O"};
};

class potentialContacts {
    /* Identifies pairs of residues that are positioned such that they could form a contact
     
     This class is designed exclusively for the scenario where we are considering a target protein and binder.
     In this case, the residue identity of the target protein residues is fixed.
     */

    public:
        potentialContacts() {};

        potentialContacts(vector<Residue*> _targetResidues, vector<Residue*> _binderResidues) {
            setTargetResidues(_targetResidues);
            setBinderResidues(_binderResidues);
        }

        void load2DProbabilityDensities(string pathToJSON);

        void setTargetResidues(vector<Residue*> _targetResidues);
        void setBinderResidues(vector<Residue*> _binderResidues);

        vector<pair<Residue*,Residue*>> getContacts(bool simpleDistanceCheck = false);

        bool isPotentialSSContact(Residue* Ri, Residue* Rj, mstreal pDensityThresh = 0.05);
        bool isPotentialSBContact(Residue* Ri, Residue* Rj, mstreal pDensityThresh = 0.05, bool checkReverseDirection = false);
        bool isPotentialBBContact(Residue* Ri, Residue* Rj);

        /**
         * @brief Constructs a idealized Cb coordinates given the backbone atoms of a residue
         * 
         * @param R A residue with all backbone heavy atoms
         * @return CartesianPoint Coordinates of idealized Cb atom
         */
        CartesianPoint getCbFromRes(Residue* R);

        mstreal getCaDistance(Residue* Ri, Residue* Rj);

        /**
         * @brief Provides a measure of how much the Ca-Cb vectors of two residues are pointed at one another
         * 
         * @param Ri 
         * @param Rj 
         * @return mstreal 
         */
        mstreal getNormalizedCbDistance(Residue* Ri, Residue* Rj);

        /**
         * @brief Provides a measure of how much the Ca-Cb vector of Ri points at the Ca atom of Rj. Note that this is a directional property. f(x,y) != f(y,x)
         * 
         * @param Ri The residue that will be used to construct a Ca-Cb vector
         * @param Rj The Ca atom of this residue will be used to create a RiCa-RjCa vector
         * @return mstreal 
         */
        mstreal getCaCbtoRiCaRjCaAngle(Residue* Ri, Residue* Rj);

    protected:
        bool ignoreContact(Residue* Ri, Residue* Rj);

    private:
        set<string> bbAtomNames = {"N","CA","C","O"};
        checkVDWRadii checkRad;

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
        mstreal defaultDistanceToCheck = 15.0;

        bool loadedPDensity = false;
        map<res_t,twoDimHistogram> sidechainSidechainContactProbabilityDensity;
        map<res_t,twoDimHistogram> sidechainBackboneContactProbabilityDensity;

        bool sameChain; // if true, will check whether nearby residues are in the same chain
        int ignoreDistance = 8; //interacting residues this close in the chain will be ignored
};

#endif