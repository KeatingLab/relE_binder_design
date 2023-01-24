#ifndef _ALIGNFRAMES_H
#define _ALIGNFRAMES_H

#include "mstrotlib.h"
#include "mstsequence.h"
#include "mstsystem.h"
#include "msttypes.h"
#include "msttransforms.h"

#include "residueframe.h"

class frameDB;
class resPair;
class resPairDB;

class augmentedStructureDB {
    public:
        augmentedStructureDB() {};
        augmentedStructureDB(string dbPath) {
            dbName = MstSys::splitPath(dbPath,1);
            readDBFile(dbPath);
        }

        ~augmentedStructureDB() {
            for (augmentedStructure *aS : targets) delete aS;
        }

        int numTargets() { return targets.size(); }
        const augmentedStructure& getTarget(int i);
        Residue* getResidue(int target_i, int res_i);
        residueFrame* getResidueFrame(int target_i, int res_i);

        void addTarget(augmentedStructure *S);
        void addTarget(const augmentedStructure &S);
        void readDBFile(string dbPath);
        void writeDBFile(string dbPath);

        const set<int>& getContacts(int target_i, int res_i, bool BBint = false);

        void setVDWContacts(int target_i, const map<int, set<int>> &_vdwContacts, bool BBint = false)
        {
            if (BBint) vdwContactsBB.insert(pair<int,map<int, set<int> > >(target_i,_vdwContacts));
            else vdwContacts.insert(pair<int,map<int, set<int> > >(target_i,_vdwContacts));
        }

    protected:

    private:
        string dbName = "";
        int version = 1;
        vector<augmentedStructure *> targets;

        map<int, map<int, set<int>>> vdwContacts; // All contact types, except backbone
        map<int, map<int, set<int>>> vdwContactsBB;
};

class mobileFrame : public residueFrame {
public:
    mobileFrame() {};
    mobileFrame(Residue* R, int _target, int _res_i, int _res_j, res_t _res_i_aa) : residueFrame(R), target(_target), res_i(_res_i), res_j(_res_j), res_i_aa(_res_i_aa)
    {
        res_j_aa = SeqTools::aaToIdx(R->getName());
    }

    void writeData(ostream& ofs);
    void readData(istream& ifs);

    int getTarget() {return target;}
    int getResI() {return res_i;}
    int getResJ() {return res_j;}
    res_t getResIIndex() {return res_i_aa;}
    res_t getResJIndex() {return res_j_aa;}

    string getName() {
        return MstUtils::toString(target) + "_" + MstUtils::toString(res_i) + "_" + SeqTools::idxToTriple(res_i_aa) + "_" + MstUtils::toString(res_j) + "_" + SeqTools::idxToTriple(res_j_aa);
    }

protected:
    void writeCartesianPointToBin(const CartesianPoint& p, ostream& ofs);
    void readCartesianPointFromBin(CartesianPoint& p, istream& ifs);

private:
    int target = -1;
    int res_i = -1;
    int res_j = -1;
    res_t res_i_aa = SeqTools::unknownIdx();
    res_t res_j_aa = SeqTools::unknownIdx();

    friend class augmentedStructureDB;
};

class alignInteractingFrames
{
public:
    alignInteractingFrames(bool _verbose = false, int _flanking_res = -1, mstreal _subsample_flanking = 0.1) : verbose(_verbose), flanking_res(_flanking_res), subsample_flanking(_subsample_flanking) {};
    
    /**
     * @brief Construct a new align Frames object by reading augmented structures from DB
     * 
     * @param dbPath path to the augmented structures DB
     */
    alignInteractingFrames(string dbPath, bool _verbose = false, int _flanking_res = -1, mstreal _subsample_flanking = 0.1) : db(dbPath), refFrame(Frame()), verbose(_verbose), flanking_res(_flanking_res), subsample_flanking(_subsample_flanking) {
        if (flanking_res >= 0 ) cout << "Warning: if flanking_res does not match the value used when defining VDW contacts in the structure database (originally 8), then some residue pairs could be double counted" << endl;
        cout << "Flanking residues: " << flanking_res << endl;
        cout << "Flanking residue subsample rate: " << subsample_flanking << endl; 
    }
    
    ~alignInteractingFrames() {
        for (mobileFrame* frame : allInteractingFrames) delete frame;
    }

    void setAA(string resName);
    void setRefFrame(const residueFrame &_rFrame) { refFrame = _rFrame; }
    void setHomThresh(mstreal _homThresh) {homologyThreshold = _homThresh;}

    void findMobileFrames();

    int getNumInteracting() { return allInteractingFrames.size(); }

    Structure *getAlignedInteractingRes(int i);
    Structure *getAlignedInteractingRes(mobileFrame* frame);
    void writeAlignedInteractingResToPDB(string pdbPath, mstreal subsampleRate = 0.01);
    void writeMobileFramesToBin(frameDB* frameBin);
    void writeResiduePairsToBin(resPairDB* rPBin);
    void writeInteractionData(string pathPrefix);

    // mobileFrame* getInteractingResidueFrame(int i);
    // void getAllInteractingResidueFrames();

    augmentedStructureDB& getDB() {return db;}

    bool isQueryHomologousToMatchInDB(Residue* query, Residue* reference, mobileFrame* mFrame);

protected:
    Structure *changeFrameToRef(Structure* interactingResiduePair, mobileFrame* frame);
    Structure *constructStructureFromResiduePair(Residue *Ri, Residue *Rj, mobileFrame* frame);

    /**
     * @brief Define a window around each residue, extract sequence, and find the identity
     * 
     * @param R1 
     * @param R2 
     * @return pair<int,int> The total sequence length and number of identical residues, respectively
     */
    pair<int,int> getSequenceIdentity(Residue* R1, Residue* R2);

private:
    augmentedStructureDB db;
    int windowSize = 30;
    mstreal homologyThreshold = 0.6;
    int flanking_res; // if greater than 0, will automatically include flanking residues in the chain
    mstreal subsample_flanking; // If greater than 0, will randomly subsample flanking residues to reduce memory footprint
    bool verbose;

    map<string, string> aaConversions;

    res_t aa = SeqTools::unknownIdx();
    string aaName = "UNK";
    Frame refFrame;

    vector<mobileFrame*> allInteractingFrames;
    vector<resPair> allInteractingRes;
    

    Transform tf;

    map<res_t,map<int,int> > interactionData;
};

class frameDB {
    public:
        frameDB(string _frameDBPath, bool _read = true) : frameDBPath(_frameDBPath), readMode(_read) {
            cout << "read mode: " << readMode << "\tfrom file: " << frameDBPath << endl;
            openFileStream();
        }

        frameDB(const frameDB& other) : frameDBPath(other.frameDBPath), readMode(other.readMode) {
            // MstUtils::assert(other.readMode, "Copying write-only frameDB file not supported");
            cout << "Opening file stream for copy of frame DB: " << frameDBPath << endl;
            openFileStream();
        }

        ~frameDB() {
            // if ((!readMode)&&(frameAdded)) MstUtils::writeBin(fs,'E');
            fs.close();
        }

        bool hasNext();
        void skip();
        void reset();

        mobileFrame* next();
        vector<mobileFrame*> loadAllFrames();

        void appendFrame(mobileFrame* rF);

    protected:
        void openFileStream();
        mobileFrame* readNextFileSection();

    private:
        string frameDBPath = "";
        int version = 1;
        bool readMode = true;
        bool frameAdded = false;

        fstream fs;
};

class resPair {
    public:
        resPair() {};
        resPair(Residue* Ri, Residue* Rj, int _target = -1);
        ~resPair() {
            if (ownsAtoms) for (Atom* A: resPairAtoms) delete A;
        }
        resPair(const resPair &rP) {
            target = rP.target;
            res_i = rP.res_i;
            res_j = rP.res_j;
            res_i_aa = rP.res_i_aa;
            res_j_aa = rP.res_j_aa;
            CaDistance = rP.CaDistance;
            resFrameBasisVectorAngles = rP.resFrameBasisVectorAngles;
            bbAtomDistances = rP.bbAtomDistances;
            ownsAtoms = rP.ownsAtoms;
            if (ownsAtoms) {
                for (Atom* A: rP.resPairAtoms) resPairAtoms.push_back(new Atom(*A));
            } else {
                resPairAtoms = rP.resPairAtoms;
            }
        }

        resPair& operator=(const resPair& rP) {
            if (this == &rP) return *this;
            target = rP.target;
            res_i = rP.res_i;
            res_j = rP.res_j;
            res_i_aa = rP.res_i_aa;
            res_j_aa = rP.res_j_aa;
            CaDistance = rP.CaDistance;
            resFrameBasisVectorAngles = rP.resFrameBasisVectorAngles;
            bbAtomDistances = rP.bbAtomDistances;
            ownsAtoms = rP.ownsAtoms;
            if (ownsAtoms) {
                for (Atom* A: rP.resPairAtoms) resPairAtoms.push_back(new Atom(*A));
            } else {
                resPairAtoms = rP.resPairAtoms;
            }
            return *this;
        }

        void writeData(ostream& ofs);
        void readData(istream& ifs);

        int getTarget() {return target;}
        int getResI() {return res_i;}
        int getResJ() {return res_j;}
        res_t getResIAAIndex() {return res_i_aa;}
        res_t getResJAAIndex() {return res_j_aa;}
        mstreal getCaDistance() {return CaDistance;}
        CartesianPoint getbbAtomDistances() {return bbAtomDistances;}
        CartesianPoint getAllbbAtomDistances();
        CartesianPoint getAngles() {return resFrameBasisVectorAngles;}
        vector<Atom*> getAtoms() {return resPairAtoms;} //Careful! The atoms are stripped of all info other than coordinates

        string getName() {
            return MstUtils::toString(target) + "_" + MstUtils::toString(res_i) + "_" + SeqTools::idxToTriple(res_i_aa) + "_" + MstUtils::toString(res_j) + "_" + SeqTools::idxToTriple(res_j_aa);
        }

protected:
    // void computeDistancesFromBBAtoms();
    void computeInternalRepresentation();

    void writeAtomToBin(Atom* A, ostream& ofs);
    Atom* readAtomFromBin(istream& ifs);

private:
    int target = -1;
    int res_i = -1;
    int res_j = -1;
    res_t res_i_aa = SeqTools::unknownIdx();
    res_t res_j_aa = SeqTools::unknownIdx();

    /*
    Ri: N, Ca, C, O; Rj: N, Ca, C, O (8 total, always in that order)
    If the resPair was directly created from existing residues, then it does not manage it's own memory (this is to avoid unecessary duplication)
    */
    vector<Atom*> resPairAtoms;
    mstreal CaDistance; // Ri:Ca - Rj:Ca
    CartesianPoint resFrameBasisVectorAngles; // the angle between each pair of residue frame basis vectors
    CartesianPoint bbAtomDistances; // Ri:N - Rj:N, Ri:Ca-Rj:Ca, Ri:C-Rj:C (3 distances. Oxygen not considered)
    bool ownsAtoms = false;

    friend class augmentedStructureDB;
};

class resPairDB {
    public:
        resPairDB(string _resPairDBPath, bool _read = true) : resPairDBPath(_resPairDBPath), readMode(_read) {
            cout << "read mode: " << readMode << "\tfrom file: " << resPairDBPath << endl;
            openFileStream();
        }

        resPairDB(const resPairDB& other) : resPairDBPath(other.resPairDBPath), readMode(other.readMode) {
            // MstUtils::assert(other.readMode, "Copying write-only frameDB file not supported");
            cout << "Opening file stream for copy of frame DB: " << resPairDBPath << endl;
            openFileStream();
        }

        ~resPairDB() {
            // if ((!readMode)&&(frameAdded)) MstUtils::writeBin(fs,'E');
            fs.close();
        }

        bool hasNext();
        void skip();
        void reset();

        resPair* next();
        vector<resPair*> loadAllResPairs();

        void appendResPair(resPair* rP);

    protected:
        void openFileStream();
        resPair* readNextFileSection();

    private:
        string resPairDBPath = "";
        int version = 1;
        bool readMode = true;
        bool resPairAdded = false;

        fstream fs;
};

#endif