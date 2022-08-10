#ifndef _SEARCHRESPAIRS_H
#define _SEARCHRESPAIRS_H

#include "msttypes.h"

#include "alignframes.h"
#include "hashframes.h"
#include "residueframe.h"

class findResPairs {
    public:
        findResPairs(string resPairDBPath, mstreal _maxDistance = 0.5, mstreal _maxRMSD = 0.5);

        ~findResPairs() {
            for (resPair* rP : allResPairs) delete rP;
        }

        void setQuery(Residue* Ri, Residue* Rj) {
            queryRP = resPair(Ri,Rj);
        }

        int searchForMatches(bool verbose = true);

        vector<resPair*> getMatches() {return verifiedMatches;}
        int getNumMatchesWithResidueType(bool Ri = true);

    protected:

    private:
        resPairDB DB;
        vector<resPair*> allResPairs;
        ProximitySearch PS;

        resPair queryRP;
        mstreal maxDistance;
        mstreal maxRMSD;

        vector<resPair*> matches;
        vector<resPair*> verifiedMatches;

        RMSDCalculator calc;
};

class findResPairsCosAngle {
    public:
        findResPairsCosAngle(string resPairDBPath, mstreal _maxDistance = 0.5, mstreal maxCosAngle = 30, mstreal _maxRMSD = 0.25);

        ~findResPairs() {
            for (resPair* rP : allResPairs) delete rP;
        }

        void setQuery(Residue* Ri, Residue* Rj) {
            queryRP = resPair(Ri,Rj);
        }

        int searchForMatches(bool verbose = true);

        vector<resPair*> getMatches() {return verifiedMatches;}
        int getNumMatchesWithResidueType(bool Ri = true);

    protected:

    private:
        resPairDB DB;
        vector<resPair*> allResPairs;
        ProximitySearch PS;

        resPair queryRP;
        mstreal maxDistance;
        mstreal maxRMSD;

        vector<resPair*> matches;
        vector<resPair*> verifiedMatches;

        RMSDCalculator calc;
};


#endif