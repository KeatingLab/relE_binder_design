#ifndef _SEARCHRESPAIRS_H
#define _SEARCHRESPAIRS_H

#include "msttypes.h"

#include "alignframes.h"
#include "hashframes.h"
#include "residueframe.h"

class distance3AngleTable {
    public:
        distance3AngleTable(mstreal _maxVal, mstreal _dCut, mstreal _angleCut) : maxVal(_maxVal), dCut(_dCut), angleCut(_angleCut) {
            nBins = ceil(2*((maxVal)/dCut));
            cout << "nBins distance: " << nBins << endl;
            binLength = maxVal/mstreal(nBins-1); // should be an extra bin at the end
            cout << "binLength: " << binLength << endl;
            for (int i = 0; i < nBins; i++) {
                // cout << "angleCut: " << angleCut << endl;
                cout << "nBins angles: " << 2*ceil(180/angleCut) << endl;
                bins.emplace_back(0,0,0,180,180,180,ceil(2*180/angleCut));
            }
        }

        void addResPairToTable(resPair* rP, int tag) {
            if (rP->getCaDistance() > 0.0) {
                // int idx = val2Idx(rP->getCaDistance());
                // cout << "bin idx: " << val2Idx(rP->getCaDistance()) << endl;
                // cout << bins[idx].getXHigh() << " " << bins[idx].getYHigh() << " " << bins[idx].getZHigh() << endl;
                bins[val2Idx(rP->getCaDistance())].addPoint(rP->getAngles(),tag);
                idx2distance[tag] = rP->getCaDistance();
            }
        }

        vector<int> searchResPair(resPair* rP);

        vector<int> getBinIdxRange(mstreal val, mstreal tol) {
            int lowerBin = val2Idx(val-tol);
            int higherBin = val2Idx(val+tol);
            cout << "lowerBin: " << lowerBin << " higherBin: " << higherBin << endl;
            vector<int> range;
            for (int i = lowerBin; i <= higherBin; i++) range.push_back(i);
            return range;
        }
    
    protected:
        int val2Idx(mstreal val) {
            if (bins.empty()) MstUtils::error("Bins are empty","lookup1D::val2Idx");
            if (val > maxVal) MstUtils::error("Provided value: "+MstUtils::toString(val)+" greater than the maximum value in the table","lookup1D::val2Idx");
            return floor(val/binLength);
        }

    private:
        mstreal maxVal;
        int nBins;
        mstreal binLength;
        vector<ProximitySearch> bins;
        map<int,mstreal> idx2distance;

        mstreal dCut;
        mstreal angleCut;
};

class findResPairs {
    public:
        findResPairs(string resPairDBPath, mstreal _dCut = 0.5, mstreal _angleCut = 30, mstreal _rmsdCut = 0.5);

        ~findResPairs() {
            for (resPair* rP : allResPairs) delete rP;
            // if (resPairMap != nullptr) delete resPairMap;
        }

        void setQuery(Residue* Ri, Residue* Rj) {
            queryRP = resPair(Ri,Rj);
        }

        resPair& getQuery() {return queryRP;}

        int searchForMatches(bool verbose = true);

        resPair* getMatch(int matchIdx);
        vector<resPair*> getMatches() {return verifiedMatches;}
        int getNumMatchesWithResidueType(bool Ri = true);
        vector<int> getNumMatchesByAAType(bool Ri = true);

    protected:

    private:
        resPairDB DB;
        vector<resPair*> allResPairs;
        
        // distance3AngleTable* resPairMap = nullptr;
        ProximitySearch PS;

        resPair queryRP;
        mstreal dCut;
        mstreal angleCut;
        mstreal rmsdCut;

        vector<resPair*> matches;
        vector<resPair*> verifiedMatches;

        RMSDCalculator calc;
};

#endif