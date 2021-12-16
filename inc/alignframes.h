#ifndef _ALIGNFRAMES_H
#define _ALIGNFRAMES_H

#include "mstsequence.h"
#include "msttypes.h"
#include "msttransforms.h"

#include "residueframe.h"

struct interactingRes
{
    interactingRes(int _target, int _res_i, int _res_j) : target(_target), res_i(_res_i), res_j(_res_j) {}
    int target;
    int res_i;
    int res_j;
};

class alignFrames
{
public:
    alignFrames(){};
    /**
     * @brief Construct a new align Frames object by reading augmented structures from DB
     * 
     * @param dbPath path to the augmented structures DB
     */
    alignFrames(string dbPath)
    {
        readDBFile(dbPath);
    }

    ~alignFrames()
    {
        for (augmentedStructure *aS : targets)
            delete aS;
    }

    void setAA(string resName)
    {
        aa = SeqTools::aaToIdx(resName);
        aaName = SeqTools::idxToTriple(aa);
    }
    void setRefFrame(const residueFrame &_rFrame) { rFrame = _rFrame; }

    void findInteractingRes();

    int getNumInteracting() { return allInteractingRes.size(); }

    Structure *getAlignedInteractingRes(int i);
    void writeAlignedInteractingResToPDB(string pdbPath);

    int numTargets() { return targets.size(); }
    void addTarget(augmentedStructure *S);
    void addTarget(const augmentedStructure &S);
    void writeDBFile(string dbPath);

    void setVDWContacts(const map<int, map<int, set<int>>> &_vdwContacts)
    {
        vdwContacts = _vdwContacts;
    }

protected:
    void readDBFile(string dbPath);

    Structure *constructStructureFromResiduePair(Residue *Ri, Residue *Rj, const interactingRes &data);

private:
    vector<augmentedStructure *> targets;
    map<int, map<int, set<int>>> vdwContacts;
    map<string, string> aaConversions;

    res_t aa = SeqTools::unknownIdx();
    string aaName = "UNK";
    residueFrame rFrame;

    vector<interactingRes> allInteractingRes;

    Transform tf;
};

#endif