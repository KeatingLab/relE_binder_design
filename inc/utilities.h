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
        static vector<Atom*> getBackboneAtoms(Chain* C) {
            vector<Atom*> allBackboneAtoms;
            for (Residue* R : C->getResidues()) {
                MstUtils::assert(RotamerLibrary::hasFullBackbone(R),"error: assertion failed. All residues in chain must have full backbones","MiscTools::getBackboneAtoms");
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

#endif