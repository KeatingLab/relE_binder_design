#ifndef _SCOREINTERFACE_H
#define _SCOREINTERFACE_H

#include "mstsequence.h"
#include "msttypes.h"

#include "alignframes.h"
#include "hashframes.h"
#include "residueframe.h"

class interfaceScorer {
    public:
        interfaceScorer(string frameDBPath, mstreal _posCut = -1.0, mstreal _oriCut = -1.0);

        ~interfaceScorer() {
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

        void scoreInterface(vector<pair<Residue*,Residue*>> _interfaceResidues = {}, alignInteractingFrames* alignF = nullptr);

        void scoreResiduePair(pair<Residue*,Residue*> resPair, alignInteractingFrames* alignF = nullptr);

        void writeContactScoresToFile(string name, bool append = true);
        void writeContactMatchesToFile(string name, bool append = true);
        void writeMatchStructures(string name, alignInteractingFrames& alignF);

        enum property {COVERAGE, SCORE};
        void writeContactPropertyToFile(string name, interfaceScorer::property toWrite);

        vector<Chain*> getBinderChains() {return binderChains;}
        vector<Chain*> getTargetChains() {return proteinChains;}

    protected:
        void setParams();

        void defineInterface();

        residueFrame* defineQuery(Residue* global, Residue* mobile);

    private:
        set<string> aaTypes;
        map<res_t,frameTable*> frameTables;
        // search parameters
        mstreal posCut = 1.0; //angstroms
        mstreal oriCut = 45.0; //degrees

        augmentedStructure *complex;
        vector<Chain*> binderChains;
        vector<Chain*> proteinChains;
        MstTimer timer;

        vector<pair<Residue*,Residue*>> interfaceResidues; // target, binder residue
        vector<set<mobileFrame*>> matchingFrames;
        vector<mstreal> searchTimes;

        fstream* contact_info_out = nullptr;
        fstream* match_info_out = nullptr;
};

#endif