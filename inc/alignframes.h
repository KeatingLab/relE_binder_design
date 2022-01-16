#ifndef _ALIGNFRAMES_H
#define _ALIGNFRAMES_H

#include "mstsequence.h"
#include "mstsystem.h"
#include "msttypes.h"
#include "msttransforms.h"

#include "residueframe.h"

class frameDB;

class proteinFrameDB {
    public:
        proteinFrameDB() {};
        proteinFrameDB(string dbPath) {
            dbName = MstSys::splitPath(dbPath,1);
            readDBFile(dbPath);
        }

        ~proteinFrameDB() {
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

        const set<int>& getContacts(int target_i, int res_i);

        void setVDWContacts(int target_i, const map<int, set<int>> &_vdwContacts)
        {
            vdwContacts.insert(pair<int,map<int, set<int> > >(target_i,_vdwContacts));
        }

    protected:

    private:
        string dbName = "";
        int version = 1;
        vector<augmentedStructure *> targets;

        map<int, map<int, set<int>>> vdwContacts;


};

class mobileFrame : public residueFrame
{
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

    friend class proteinFrameDB;
};

class alignInteractingFrames
{
public:
    alignInteractingFrames() {};
    
    /**
     * @brief Construct a new align Frames object by reading augmented structures from DB
     * 
     * @param dbPath path to the augmented structures DB
     */
    alignInteractingFrames(string dbPath) : db(dbPath), refFrame(Frame()) {}

    void setAA(string resName);
    void setRefFrame(const residueFrame &_rFrame) { refFrame = _rFrame; }
    void setHomThresh(mstreal _homThresh) {homologyThreshold = _homThresh;}

    void findMobileFrames();

    int getNumInteracting() { return allInteractingFrames.size(); }

    Structure *getAlignedInteractingRes(int i);
    Structure *getAlignedInteractingRes(mobileFrame* frame);
    void writeAlignedInteractingResToPDB(string pdbPath, mstreal subsampleRate = 0.01);
    void writeMobileFramesToBin(frameDB* frameBin);
    void writeInteractionData(string pathPrefix);

    // mobileFrame* getInteractingResidueFrame(int i);
    // void getAllInteractingResidueFrames();

    proteinFrameDB& getDB() {return db;}

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
    proteinFrameDB db;
    int windowSize = 30;
    mstreal homologyThreshold = 0.6;

    map<string, string> aaConversions;

    res_t aa = SeqTools::unknownIdx();
    string aaName = "UNK";
    Frame refFrame;

    vector<mobileFrame*> allInteractingFrames;

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
            MstUtils::assert(other.readMode, "Copying write-only frameDB file not supported");
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

#endif