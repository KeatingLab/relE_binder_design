#ifndef _GENERATESEEDS_H
#define _GENERATESEEDS_H

#include "mstfasst.h"
#include "mstmagic.h"
#include "mstrotlib.h"
#include "mstsystem.h"
#include "msttypes.h"

#include "residuecontact.h"
#include "freesasaext.h"
#include "utilities.h"
#include "utilitiesio.h"

using namespace std;


struct seedGenParams {
    int targetFlankRes = 1;
    int maxNumMatches = 1000;
    mstreal RMSDCutoff = 0.5;
    bool seqConst = true;
    int seedFlankRes = 2;
    bool writeToPDB = false;
    string contactData = "";
    mstreal minContsPerRes = 4.0; // potential contacts per seed res
};

class seedGenerator {
    public:
        seedGenerator(const Structure& target, string fasstDBPath, seedGenParams params = seedGenParams());
 
        ~seedGenerator() {
            for (auto frag : allFragmentsData) delete frag.fragment;
        }
        
        void defineBindingSiteAsSurfaceRes(mstreal rSASAthresh = 0.1);
        void setBindingSite(vector<Residue*> bindingSiteRes) {
            bindingSite = bindingSiteRes;
            std::cout << "Set the binding site to " << bindingSite.size() << " residues" << endl;
        }

        string generateSeeds();

        string getTargetName() {return targetName;}

        void writeBindingSiteFragments(string path);

        struct fragmentData {
            public:
                seedGenerator* parent = nullptr;
                int cenResIdxInFrag = -1;
                Residue* fragCenResInParent = nullptr;
                Structure* fragment = nullptr;
                fasstSolutionSet matches;

                string getName() {
                    if (name == "") {
                        // parentStructure_cenResChainIDNum
                        name = parent->getTargetName() + "_" + fragCenResInParent->getChainID() 
                        + MstUtils::toString(fragCenResInParent->getNum());
                    }
                    return name;
                }

            private:
                string name = "";
        };

        string getSeedName(int targetID, int ntermResIdx, int resLen) {
            string seedName = MstUtils::toString(targetID) + "-" + MstUtils::toString(ntermResIdx) + "-" + MstUtils::toString(resLen);
            return seedName;
        }

    protected:

        void getBindingSiteFragments();

        void findStructuralMatches();

        void generateSeedsFromMatches(seedBinaryFile& seedBin, fstream& fragOut, fstream& seedOut, multiPDBFile* pdbOut = nullptr);

        bool seedTargetClash(Structure* seed);

    private:
        const Structure& target;
        string targetName = "";
        vector<Atom*> targetBackboneAtoms = {}; // All atoms of core residues and all backbone atoms of surface residues
        ProximitySearch targetAtomsPS;
        checkVDWRadii checker;

        potentialContacts potConts;

        vector<Residue*> bindingSite = {};

        FASST F;

        seedGenParams params;
        vector<fragmentData> allFragmentsData;

        MstTimer timer;
};

#endif