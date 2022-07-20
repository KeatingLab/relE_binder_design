#ifndef _SCOREINTERFACE_H
#define _SCOREINTERFACE_H

#include "mstsequence.h"
#include "msttypes.h"

#include "alignframes.h"
#include "fragmentdb.h"
#include "hashframes.h"
#include "residueframe.h"
#include "utilities.h"

class interfaceSearch {
    public:
        interfaceSearch(string frameDBPath, mstreal _posCut = -1.0, mstreal _oriCut = -1.0);

        ~interfaceSearch() {
            contact_info_out->close();
            if (contact_info_out != nullptr) delete contact_info_out;
            match_info_out->close();
            if (match_info_out != nullptr) delete match_info_out;

            for (auto const it : frameTables) delete it.second;
        }

        void loadStructure(augmentedStructure *_complex, vector<string> targetChainIDs, vector<string> binderChainIDs) {
            complex = _complex;

            for (string chainID : targetChainIDs) {
                Chain* C = complex->getChainByID(chainID);
                if (C == NULL) MstUtils::error("No chain with ID: "+chainID+ " found in complex: "+complex->getName(),"interfaceScorer::loadStructure");
                proteinChains.push_back(C);
            }
            for (string chainID : binderChainIDs) {
                Chain* C = complex->getChainByID(chainID);
                if (C == NULL) MstUtils::error("No chain with ID: "+chainID+ " found in complex: "+complex->getName(),"interfaceScorer::loadStructure");
                binderChains.push_back(C);
            }
        }

        void setInterface(vector<pair<Residue*,Residue*>> _interfaceResidues) {
            interfaceResidues = _interfaceResidues;
            for (pair<Residue*,Residue*> contact : interfaceResidues) {
                Residue* targetRes = contact.first;
                targetResidues.insert(targetRes);
            }
            cout << interfaceResidues.size() << " interface contacts" << endl;
            cout << targetResidues.size() << " target residues" << endl;
        }

        residueFrame* defineQuery(Residue* global, Residue* mobile) {
            residueFrame* refFrame = complex->getResidueFrame(global->getResidueIndex());
            residueFrame* mobileFrame = complex->getResidueFrame(mobile->getResidueIndex());
            return mobileFrame->frameRelativeToOther(*refFrame);
        }

        vector<Chain*> getBinderChains() {return binderChains;}
        vector<Chain*> getTargetChains() {return proteinChains;}

        void searchInterface(alignInteractingFrames* alignF = nullptr);

        void searchResiduePair(pair<Residue*,Residue*> resPair, alignInteractingFrames* alignF = nullptr);

        void writeContactScoresToFile(string name, bool append = true);
        void writeContactMatchesToFile(string name, bool append = true);
        void writeMatchStructures(string name, alignInteractingFrames& alignF);

        void writeContactPropertyToFile(string name);

    protected:
        void defineInterface();

    private:
        set<string> aaTypes;

        augmentedStructure *complex;
        vector<Chain*> binderChains;
        vector<Chain*> proteinChains;

        vector<pair<Residue*,Residue*>> interfaceResidues; // target, binder residue
        set<Residue*> targetResidues;

        mstreal posCut = 0.25; //angstroms
        mstreal oriCut = 15.0; //degrees

        map<res_t,frameTable*> frameTables;

        MstTimer timer;

        vector<set<mobileFrame*>> matchingFrames;
        vector<mstreal> searchTimes;
 
        fstream* contact_info_out = nullptr;
        fstream* match_info_out = nullptr;
};

struct binderScorerParams {
    string frameDBPath = ""; //path to database of mobile frames aligned to a global frame
    string potentialContactsJSONPath = ""; //path to JSON file with 2D probability density
    mstreal posCut = 0; //angstroms
    mstreal oriCut = 0; //degrees
    bool renormalizeProbabilities = true;
    bool verbose = false;
};

class binderScorer {
    public:
        binderScorer(const binderScorerParams& params, augmentedStructure& _target);

        binderScorer(const binderScorerParams& params, augmentedStructure& complex, string binderChainIDsString, string targetChainIDsString);

        ~binderScorer() {
            if (complexMode) delete binder;
            contact_info_out->close();
            if (binder_info_out != nullptr) delete binder_info_out;
            if (contact_info_out != nullptr) delete contact_info_out;
            if (clash_info_out != nullptr) delete clash_info_out;
            for (auto const it : frameTables) delete it.second;
        }

        void setBinder(augmentedStructure* _binder);
        augmentedStructure* getBinder() {return binder;}
        void setTargetBindingSiteResidues(vector<Residue*> sel);
        void defineTargetBindingSiteResiduesByrSASA(mstreal relSASAthreshold = 0.05);

        int defineInterfaceByPotentialContacts();
        // void manuallySetInterface(vector<pair<Residue*,Residue*>> contacts);

        vector<Residue*> getTargetResidues() {return target.getResidues();}
        vector<Residue*> getBinderResidues() {return binder->getResidues();}

        residueFrame* defineQuery(Residue* global, Residue* mobile) {
            residueFrame* refFrame = target.getResidueFrame(global->getResidueIndex());
            residueFrame* mobileFrame = binder->getResidueFrame(mobile->getResidueIndex());
            return mobileFrame->frameRelativeToOther(*refFrame);
        }

        mstreal scoreInterface();

        int countDesignableContacts() {return pConts.getNumContacts();}
        int countNonDesignableContacts() {return pConts.getNumNonDesignablePairs();}

        void writeBinderScoresToFile(bool append = true);
        void writeContactScoresToFile(bool append = true);
        void writeContactPropertyToFile(string dirPath = "");
        void writeResidueClashesToFile(bool append = true);
        void writeResiduesClashesToSimpleFile(string dirPath = "");
        void writePSSM(string dirPath = "");

    protected:
        void prepareVoxelGrids(string frameDBPath);
        void defineInterfaceUsingVDWContacts();
        // void setFrameProbabilityTables();

        int residueFrameCounts(pair<Residue*,Residue*> resPair);

        mstreal logProb(mstreal numerator, mstreal denominator, mstreal pseudocount = 1.0);
        mstreal logCounts(int counts, mstreal pseudocount = 1.0) {return -log(counts+pseudocount);}

        vector<mstreal> getAADistAtPos(Residue* binderRes);

        /**
         * @brief Compute the joint probability from multiple amino acid distributions
         * 
         * @param aaDists The inner vector is the probability distribution, the outer vector is each distribution
         * @return vector<mstreal> 
         */
        vector<mstreal> getJointProbability(vector<vector<mstreal>> aaDists);

    private:
        string targetName = "";
        augmentedStructure target;
        Structure targetBackbone;
        augmentedStructure* binder = nullptr;
        set<string> aaTypes;
        bool verbose = false;

        vector<Chain*> binderChains;
        vector<Chain*> proteinChains;

        potentialContacts pConts;
        vector<pair<Residue*,Residue*>> interfaceResidues; // target, binder residue
        vector<pair<Residue*,Residue*>> nonDesignableInterfaceResidues; // target, binder
        set<Residue*> targetResidues; //subset of residues in the target found in the binding site

        mstreal nonDesignableContactPenalty = 10.0;
        mstreal designabilityScore = 0.0;
        int numNonDesignable = 0;

        mstreal posCut = 1.0; //angstroms
        mstreal oriCut = 10.0; //degrees

        // bool renormalizeProbabilities = true;
        map<res_t,frameTable*> frameTables;

        MstTimer timer;

        // vector<mstreal> probabilities;
        // vector<pair<int,int>> countAndNorm;
        vector<int> interfaceCounts;
        vector<mstreal> interfaceScores;
        vector<mstreal> searchTimes;

        fstream* binder_info_out = nullptr;
        fstream* contact_info_out = nullptr;
        fstream* clash_info_out = nullptr;
        // fstream* match_info_out = nullptr;

        bool complexMode = false; //indicates that we're scoring a real complex
};

class binderBackboneScorer {
    public:
        binderBackboneScorer(string segmentGraphPath, string _name) : segmentSearcher(segmentGraphPath), name(_name) {;}
        ~binderBackboneScorer() {if (info_out.is_open()) info_out.close();}

        mstreal scoreBackbone(Structure* S);

    protected:
        void openFile();

    private:
        searchSegments segmentSearcher;

        string name;
        fstream info_out;
};

#endif