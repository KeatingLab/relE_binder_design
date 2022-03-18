#ifndef _HASHFRAMES_H
#define _HASHFRAMES_H

#include <math.h>

#include "mstsequence.h"
#include "msttypes.h"
#include "msttransforms.h"

#include "alignframes.h"
#include "residuecontact.h"
#include "residueframe.h"

class boundingBox {
    public:
        boundingBox(mstreal _pad = 0.0) : pad(_pad) {};

        void update(vector<mobileFrame*> frames);
        void update(Frame* frame);
        void update(const CartesianPoint& point);

        bool isWithinBounds(Frame* frame, mstreal tolerance = 0.0) const;
        bool isWithinBounds(Atom* A, mstreal tolerance = 0.0) const;
        bool isWithinBounds(const CartesianPoint& point, mstreal tolerance = 0.0) const;

        void printBounds();

        mstreal getXMin() const {return xMin;}
        mstreal getXMax() const {return xMax;}
        mstreal getYMin() const {return yMin;}
        mstreal getYMax() const {return yMax;}
        mstreal getZMin() const {return zMin;}
        mstreal getZMax() const {return zMax;}

        mstreal getXWidth() const {return xMax-xMin;}
        mstreal getYWidth() const {return yMax-yMin;}
        mstreal getZWidth() const {return zMax-zMin;}
    private:
        mstreal pad;

        mstreal xMin = std::numeric_limits<double>::max();
        mstreal xMax = std::numeric_limits<double>::min();
        mstreal yMin = std::numeric_limits<double>::max();
        mstreal yMax = std::numeric_limits<double>::min();
        mstreal zMin = std::numeric_limits<double>::max();
        mstreal zMax = std::numeric_limits<double>::min();
};

class positionHasher {
    public:
        positionHasher(boundingBox _bbox, mstreal _increment);

        int hashFrame(Frame* frame, bool strict = false);
        int hashPoint(const CartesianPoint& point, bool strict = false);

        vector<int> region(Frame* frame, mstreal cutoff);
        vector<int> region(const CartesianPoint& point, mstreal cutoff);

        int getNumBins() {return numXBins * numYBins * numZBins;}

        CartesianPoint getCenterCoordinatesOfBin(int hash);

        const boundingBox& getBoundingBox() {return bbox;}
    private:
        int hashCoordinates(const int& xBin, const int& yBin, const int& zBin);
        void getCoordinatesFromHash(int hash, int& xBin, int& yBin, int& zBin);

        boundingBox bbox;
        mstreal increment;
        int numXBins;
        int numYBins;
        int numZBins;
};

class anglesFromFrame {
    public:
        anglesFromFrame(Frame* frame);
        
        mstreal getAlpha() const {return alphaAngle;}
        mstreal getBeta() const {return betaAngle;}
        mstreal getGamma() const {return gammaAngle;}

        bool L1Check(const anglesFromFrame& other, mstreal cutoff);

    private:
        // all angles are in radians
        mstreal alphaAngle = -1.0;
        mstreal betaAngle = -1.0;
        mstreal gammaAngle = -1.0;
};

class orientationHasher {
    public:
        orientationHasher(int _numBins);

        int hashFrame(Frame* frame);
        vector<int> region(Frame* frame, mstreal cutoff);
        int getNumBinsTotal() {return numBins * numBins * numBins;}

    protected:     
        int hashAngles(const anglesFromFrame& angles);
        void getAngleBins(const anglesFromFrame& angles, int& alphaBin, int& betaBin, int& gammaBin, bool strict = true);
        int hashAngleBins(const int& alphaBin, const int& betaBin, const int& gammaBin);
        void getAngleBinsRange(const anglesFromFrame& angles, mstreal cutoff, vector<int>& alphaBins, vector<int>& betaBins, vector<int>& gammaBins);


    private:
        mstreal increment; //radians
        int numBins; //number of bins along each axis of rotation
};

class frameTable {
    public:
        frameTable(boundingBox _bbox, mstreal distIncrement, mstreal numRotBins) : posHasher(_bbox,distIncrement), rotHasher(numRotBins) {
            posLookupTable.resize(posHasher.getNumBins());
            rotLookupTable.resize(rotHasher.getNumBinsTotal());
        }

        ~frameTable() {
            for (mobileFrame* frame : allFrames) delete frame;
        }

        void insertFrame(mobileFrame* frame);

        /**
         * @brief finds all frames satisfying search criteria w.r.t. query
         * 
         * @param frame the query
         * @param distCut the distance cutoff in angstroms
         * @param angCut the orientation cutoff in degrees
         * @return int 
         */
        int countSimilarFrames(Frame* frame, mstreal distCut, mstreal angCut, vector<bool>* posHashToIgnore = nullptr);
        set<mobileFrame*> findSimilarFrames(Frame* frame, mstreal distCut, mstreal angCut, vector<bool>* posHashToIgnore = nullptr);

        void static printFrameInfo(Frame* frame);
        const boundingBox& getBoundingBox() {return posHasher.getBoundingBox();}

        int getNumFramesInPosBin(int hash) {return posLookupTable[hash].size();}
        int getTotalNumFramesInTable() {return allFrames.size();}

    protected:
        set<mobileFrame*> verify(Frame* queryFrame, const vector<mobileFrame*>& distNeighbors, mstreal distCut, const vector<mobileFrame*>& rotNeighbors, mstreal angCut);

        positionHasher posHasher;
        orientationHasher rotHasher;
    private:

        vector<mobileFrame*> allFrames;
        vector<vector<mobileFrame*> > posLookupTable;
        vector<vector<mobileFrame*> > rotLookupTable;
};

class frameProbability : public frameTable {
    public:
        frameProbability(boundingBox _bbox, mstreal distIncrement, mstreal numRotBins) : frameTable(_bbox, distIncrement, numRotBins) {}

        void setSearchParams(mstreal _distCut, mstreal _oriCut) {
            distCut = _distCut;
            oriCut = _oriCut;
        }

        bool isTargetResidueDefined(int resIdx) {return targetResDefined.find(resIdx) != targetResDefined.end();}
        void setTargetResidue(int resIdx, const Structure& target, bool clashCheck = true);

        pair<int,int> findInteractingFrameProbability(Frame* frame, int targetResIdx);
        vector<mstreal> findBinderResAADist(Frame* frame, int targetResIdx, mstreal pseudocount = 1.0);

    protected:
        vector<mstreal> aaDistFromMobileFrames(set<mobileFrame*> frames, mstreal pseudocount);

    private:
        map<int,vector<bool> > occupiedVolMap; // map[resIdx][voxelIdx]
        map<int,int> normalizingConstant;  //map[resIdx]
        set<int> targetResDefined;

        mstreal distCut = 0.25;
        mstreal oriCut = 15.0;

        mstreal clashTol = 0.95;

        checkVDWRadii atomVDWcheck;
};

#endif