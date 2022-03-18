#ifndef _UTILITIES_H
#define _UTILITIES_H

#include "mstrotlib.h"
#include "mstsequence.h"

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

#endif