#ifndef _HASHFRAMES_H
#define _HASHFRAMES_H

#include <math.h>

#include "mstsequence.h"
#include "msttypes.h"
#include "msttransforms.h"

#include "alignframes.h"
#include "residueframe.h"

class positionHasher {
    public:
        positionHasher(vector<mstreal> _bbox, mstreal _increment);

        int hashFrame(Frame* frame);

        vector<int> region(Frame* frame, mstreal cutoff);

        static void updateBoundingBox(vector<mobileFrame*> frames, vector<mstreal>& bbox, bool verbose = false);
    private:
        int hashCoordinates(const int& xBin, const int& yBin, const int& zBin);

        vector<mstreal> bbox;
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

    protected:     
        int hashAngles(const anglesFromFrame& angles);
        void getAngleBins(const anglesFromFrame& angles, int& alphaBin, int& betaBin, int& gammaBin);
        int hashAngleBins(const int& alphaBin, const int& betaBin, const int& gammaBin);
        void getAngleBinsRange(const anglesFromFrame& angles, mstreal cutoff, vector<int>& alphaBins, vector<int>& betaBins, vector<int>& gammaBins);


    private:
        mstreal increment; //radians
        int numBins; //number of bins along each axis of rotation
};

class frameTable {
    public:
        frameTable(vector<mstreal> _bbox, mstreal distIncrement, mstreal numRotBins) : posHasher(_bbox,distIncrement), rotHasher(numRotBins) {}

        void insertFrame(mobileFrame* frame);

        int countSimilarFrames(Frame* frame, mstreal distCut, mstreal angCut);
        set<mobileFrame*> findSimilarFrames(Frame* frame, mstreal distCut, mstreal angCut);

    protected:
        set<mobileFrame*> verify(Frame* queryFrame, const vector<mobileFrame*>& distNeighbors, mstreal distCut, const vector<mobileFrame*>& rotNeighbors, mstreal angCut);

    private:
        positionHasher posHasher;
        orientationHasher rotHasher;

        unordered_map<int,vector<mobileFrame*> > posMap;
        unordered_map<int,vector<mobileFrame*> > rotMap;
};

#endif