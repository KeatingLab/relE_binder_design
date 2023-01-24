#ifndef _UTILITIES_H
#define _UTILITIES_H

#include "mstrotlib.h"
#include "mstsequence.h"

class fisherYatesShuffle {
    // https://en.wikipedia.org/wiki/Fisher%E2%80%93Yates_shuffle
    // Durstenfield's modern version for sampling without replacement in linear time
    public:
        /**
         * @brief Sample values in the range [minValue,maxValue]
         * 
         * @param minValue 
         * @param maxValue 
         */
        fisherYatesShuffle(int minValue, int maxValue) {
            values = vector<int>(maxValue-minValue+1);
            iota(values.begin(),values.end(),minValue);
            numRemainingToBeSampled = maxValue-minValue+1;
        }
        fisherYatesShuffle() {}

        int numRemaining() {return numRemainingToBeSampled;}

        int sample() {
            if (values.empty()) MstUtils::error("No values to sample from","fisherYatesShuffle::sample");
            // select random position in vector and get value
            int i = MstUtils::randInt(numRemainingToBeSampled);
            int val = values[i];

            // swap in whatever value is at the end of the vector
            values[i] = values[numRemainingToBeSampled-1];

            // decrease the range of the vector from which values can be sampled
            numRemainingToBeSampled--;
            return val;
        }

    private:
        vector<int> values;
        int numRemainingToBeSampled;
};

class SeqToolsExtension {
    public:
        static bool initConstants();

        static set<string> getAANames();

        static string findEquivalentResidueInAlphabet(string aa3, bool strict = true);

        static bool AAinSet(string aa3);

        static int numAA();

    private:
        static set<string> aaNames;
        static map<string,string> equivalentAA;
};

class MiscTools {
    public:
        static vector<Atom*> getBackboneAtoms(Structure& S) {
            vector<Chain*> Cvec;
            for (int i = 0; i < S.chainSize(); i++) {
                Cvec.push_back(&S.getChain(i));
            }
            return getBackboneAtoms(Cvec);
        }

        static vector<Atom*> getBackboneAtoms(vector<Chain*> Cvec) {
            vector<Atom*> allBackboneAtoms;
            for (Chain* C : Cvec) {
                vector<Atom*> chainBackboneAtoms = getBackboneAtoms(C);
                allBackboneAtoms.insert(allBackboneAtoms.end(),chainBackboneAtoms.begin(),chainBackboneAtoms.end());
            }
            return allBackboneAtoms;
        }

        static vector<Atom*> getBackboneAtoms(Chain* C) {
            vector<Atom*> allBackboneAtoms;
            for (Residue* R : C->getResidues()) {
                if (!RotamerLibrary::hasFullBackbone(R)) MstUtils::error("All residues in chain must have full backbones","MiscTools::getBackboneAtoms");
                vector<Atom*> backboneAtoms = RotamerLibrary::getBackbone(R);
                allBackboneAtoms.insert(allBackboneAtoms.end(),backboneAtoms.begin(),backboneAtoms.end());
            }
            return allBackboneAtoms;
        }

        static vector<Atom*> getBackboneAtoms(vector<Residue*> residues) {
            vector<Atom*> allBackboneAtoms;
            for (Residue* R : residues) {
                if (!RotamerLibrary::hasFullBackbone(R)) MstUtils::error("All residues in chain must have full backbones","MiscTools::getBackboneAtoms");
                vector<Atom*> backboneAtoms = RotamerLibrary::getBackbone(R);
                allBackboneAtoms.insert(allBackboneAtoms.end(),backboneAtoms.begin(),backboneAtoms.end());
            }
            return allBackboneAtoms;
        }

        static vector<Atom*> getBackboneCaAtoms(vector<Residue*> residues) {
            vector<Atom*> backboneCaAtoms;
            for (Residue* R : residues) {
                if (!RotamerLibrary::hasFullBackbone(R)) MstUtils::error("All residues in chain must have full backbones","MiscTools::getBackboneAtoms");
                vector<Atom*> backboneAtoms = RotamerLibrary::getBackbone(R);
                for (Atom* A : backboneAtoms) if (RotamerLibrary::backboneAtomType(A) == 1) backboneCaAtoms.push_back(A); // CA defined as 1 in mstrotlib.h
            }
            return backboneCaAtoms;
        }

        struct alignment {
            mstreal rmsd = std::numeric_limits<double>::max();
            int CiResIdx = -1;
            int CjResIdx = -1;
            int length = -1;
        };

        /**
         * @brief Finds the windows of two fixed-in-space structures with the lowest backbone RMSD
         * 
         * @param Ci The first chain
         * @param Cj The second chain
         * @param kmerL The length of the window of comparison
         * @return mstreal The RMSD of the window
         */
        static alignment bestRMSD(Chain* Ci, Chain* Cj, int kmerL = 5);

        static void extractBackboneFromStructure(const Structure& source, Structure& destination);
};

class twoDimHistogram {
    public:
        twoDimHistogram() {};

        twoDimHistogram(int _dim1NBins, int _dim2NBins, mstreal _dim1Min, mstreal _dim1Max, mstreal _dim2Min, mstreal _dim2Max) : dim1NBins(_dim1NBins), dim2NBins(_dim2NBins), dim1Min(_dim1Min), dim1Max(_dim1Max), dim2Min(_dim2Min), dim2Max(_dim2Max) {
            dim1BinSize = (dim1Max-dim1Min)/mstreal(dim1NBins);
            dim2BinSize = (dim2Max-dim2Min)/mstreal(dim2NBins);
            normalizingConstant = 0;
            counts.resize(dim1NBins);
            for (int i = 0; i < dim1NBins; i++) counts[i].resize(dim2NBins,0);
        };

        void setCounts(vector<vector<int>> _counts);

        int getCounts(mstreal val1, mstreal val2);
        mstreal getDensity(mstreal val1, mstreal val2);
        mstreal getdim1Min() {return dim1Min;}
        mstreal getdim1Max() {return dim1Max;}
        mstreal getdim2Min() {return dim2Min;}
        mstreal getdim2Max() {return dim2Max;}

        void reportHistogram();

    protected:
        void getIdx(mstreal val1, mstreal val2, int& idx1, int& idx2);

    private:
        int dim1NBins;
        int dim2NBins;
        mstreal dim1Min, dim1Max;
        mstreal dim2Min, dim2Max;
        mstreal dim1BinSize;
        mstreal dim2BinSize;
        vector<vector<int>> counts;
        int normalizingConstant;
};

class lineSDF {
    // Derivation from an Inigo Quilez video https://youtu.be/PMltMdi1Wzg
    public:
        lineSDF(CartesianPoint _A, CartesianPoint _B, bool _roundedEdges = false) : A(_A), B(_B), roundedEdges(_roundedEdges) {;}

        mstreal distance(CartesianPoint P) {
            mstreal h = normalizedDistanceAlongLine(P);
            if (roundedEdges) {
                return (P - A - (B-A)*min(1.0,max(0.0,h))).norm();
            } else {
                if (h > 1.0) return std::numeric_limits<double>::max();
                if (h < 0.0) return std::numeric_limits<double>::max();
                return (P - A - (B-A)*h).norm();
            }
        }

    protected:
        mstreal normalizedDistanceAlongLine(CartesianPoint P) {
            mstreal normOfAtoB = (B-A).norm();
            return (P - A).dot(B-A)/(normOfAtoB*normOfAtoB);
        }

    private:
        CartesianPoint A;
        CartesianPoint B;
        bool roundedEdges = false; // If true, then report distance of points that are closer to the endpoints
};

#endif