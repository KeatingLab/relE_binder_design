#ifndef _SCOREINTERFACE_H
#define _SCOREINTERFACE_H

#include "mstsequence.h"
#include "msttypes.h"

#include "alignframes.h"
#include "fragmentdb.h"
#include "hashframes.h"
#include "residueframe.h"
#include "searchrespairs.h"
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
    // Paths to databases
    string frameDBPath = ""; //path to database of mobile frames aligned to a global frame
    string resPairDBPath = ""; //path to the database of interacting residue backbone atoms
    string potentialContactsJSONPath = ""; //path to JSON file with 2D probability density

    // Search parameters
    mstreal posCut = 0; //angstroms (residue frame)
    mstreal oriCut = 0; //degrees (residue frame)
    mstreal dCut = 0; //angstroms (residue backbone atoms)
    mstreal angleCut = 0; // degrees (between residue frame basis vectors)
    mstreal RMSDCut = 0; //angstroms (residue backbone atoms)

    // Scoring parameters
    mstreal nonDesignableContactPenalty = 10.0; // penalty for having a contact geometry not observed in the PDB
    mstreal noMatchesPenalty = 10.0; // penalty for having a residue backbone geometry not obsered in the PDB
    
    // Options for scoring
    bool renormalizeProbabilities = true;
    bool verbose = false;
};

// A base class that handles shared features of the binder scoring classes
class binderScorer {
    public:
        binderScorer(const binderScorerParams& _params) : params(_params) {
            aaTypes = SeqToolsExtension::getAANames();
            pConts.load2DProbabilityDensities(params.potentialContactsJSONPath);
            setBackgroundSurfaceProbabilities();
        }

        /**
         * @brief Constructor using only a target structure
         * 
         * @param params 
         * @param _target 
         */
        binderScorer(const binderScorerParams& _params, Structure& _target);
        
        /**
         * @brief Constructor that takes a structure of the complex and converts into a proper format for scoring
         * 
         * @param params 
         * @param complex 
         * @param binderChainIDsString 
         * @param targetChainIDsString 
         */
        binderScorer(const binderScorerParams& _params, Structure& complex, string binderChainIDsString, string targetChainIDsString);

        ~binderScorer() {
            if (complexMode) delete binder;
        }

        void setComplex(Structure& complex, string binderChainIDsString, string targetChainIDsString);
        void setBinder(Structure* _binder);
        void setTarget(Structure& target);
        void setTargetBindingSiteResidues(vector<Residue*> sel);
        void defineTargetBindingSiteResiduesByrSASA(mstreal relSASAthreshold = 0.05);
        
        int defineInterfaceByPotentialContacts();
        int countDesignableContacts() {return pConts.getNumContacts();}

        Structure* getBinder() {return binder;}
        vector<Residue*> getTargetResidues() {return target.getResidues();}
        vector<Residue*> getBinderResidues() {return binder->getResidues();}

    protected:
        void defineInterfaceUsingVDWContacts();
        void setBackgroundSurfaceProbabilities();

        string targetName = "";
        Structure target;
        Structure* binder = nullptr;
        set<string> aaTypes;
        map<res_t,mstreal> backgroundSurfaceProbabilities;

        vector<Chain*> binderChains;
        vector<Chain*> proteinChains;

        potentialContacts pConts;
        vector<pair<Residue*,Residue*>> interfaceResidues; // target, binder residue
        vector<pair<Residue*,Residue*>> nonDesignableInterfaceResidues; // target, binder residue
        set<Residue*> targetBindingResidues; //subset of residues in the target found in the binding site
        
        binderScorerParams params;
        MstTimer timer;
        bool complexMode = false; //indicates that a complex structure (target + binder) was provided
};

/* --- --- --- --- --- residueBackboneBinderScorer --- --- --- --- --- */

class residueBackboneBinderScorer: public binderScorer {
    public:
        residueBackboneBinderScorer(const binderScorerParams& params) : binderScorer(params), resPairSearcher(params.resPairDBPath,params.dCut,params.angleCut,params.RMSDCut) {};

        residueBackboneBinderScorer(const binderScorerParams& params, Structure& target) : binderScorer(params,target), resPairSearcher(params.resPairDBPath,params.dCut,params.angleCut,params.RMSDCut) {};

        residueBackboneBinderScorer(const binderScorerParams& params, Structure& complex, string binderChainIDsString, string targetChainIDsString) : binderScorer(params,complex,binderChainIDsString,targetChainIDsString), resPairSearcher(params.resPairDBPath,params.dCut,params.angleCut,params.RMSDCut) {};

        mstreal scoreBinder();

        void writeBinderScoresToFile(bool append = true);
        void writeContactScoresToFile(bool append = true);
        void writeTrainingDataToFile(bool append = true);

        // void writeBinderPSSMToFile(bool append = true);

    protected:
        void scoreContact(pair<Residue*,Residue*> contactingRes, int& totalMatches, int& nativeMatches, mstreal& contactScore);

    private:
        findResPairs resPairSearcher;

        mstreal binderScore = 0.0;
        vector<mstreal> binderScorePerContact;
        vector<int> numMatchesPerContact;
        vector<int> numNativeMatchesPerContact;
        vector<vector<int>> nMatchesAAPerContact;
        vector<mstreal> searchTimePerContact;
        vector<CartesianPoint> backboneAtomDistancesPerContact;

        fstream* binder_info_out = nullptr;
        fstream* contact_info_out = nullptr;
        fstream* training_info_out = nullptr;

};

// class residueFrameBinderScorer : public binderScorer {
//     public:
//         residueFrameBinderScorer(const binderScorerParams& params, augmentedStructure& _target);

//         residueFrameBinderScorer(const binderScorerParams& params, augmentedStructure& complex, string binderChainIDsString, string targetChainIDsString);

//         ~residueFrameBinderScorer() {
//             if (complexMode) delete augmentedBinder;
//             contact_info_out->close();
//             if (binder_info_out != nullptr) delete binder_info_out;
//             if (contact_info_out != nullptr) delete contact_info_out;
//             if (clash_info_out != nullptr) delete clash_info_out;
//             for (auto const it : frameTables) delete it.second;
//         }

//         residueFrame* defineQuery(Residue* global, Residue* mobile) {
//             residueFrame* refFrame = target.getResidueFrame(global->getResidueIndex());
//             residueFrame* mobileFrame = binder->getResidueFrame(mobile->getResidueIndex());
//             return mobileFrame->frameRelativeToOther(*refFrame);
//         }

//         mstreal scoreInterfaceDesignability();
//         mstreal scoreInterfaceSequenceCompatibility();

//         int countDesignableContacts() {return pConts.getNumContacts();}
//         int countNonDesignableContacts() {return pConts.getNumNonDesignablePairs();}

//         void writeBinderDesignabilityScoresToFile(bool append = true);
//         void writeContactDesignabilityScoresToFile(bool append = true);

//         void writeBinderSequenceCompatibilityScoresToFile(bool append = true);
//         void writeContactSequenceCompatibilityScoresToFile(bool append = true);

//         void writeContactPropertyToFile(string dirPath = "", bool designabilityScore = true);

//         void writeResidueClashesToFile(bool append = true);
//         void writeResiduesClashesToSimpleFile(string dirPath = "");
//         void writePSSM(string dirPath = "");

//     protected:
//         void prepareVoxelGrids(string frameDBPath);
//         // void setFrameProbabilityTables();

//         int residueFrameCounts(pair<Residue*,Residue*> resPair);

//         void residueFrameSequenceCompatibility(pair<Residue*,Residue*> resPair, int &totalMatches, int &nativeAAMatches, mstreal &sequenceCompatibilityScore);

//         mstreal logProb(mstreal numerator, mstreal denominator, mstreal pseudocount = 1.0);
//         mstreal logCounts(int counts, mstreal pseudocount = 1.0) {return -log(counts+pseudocount);}

//         vector<mstreal> getAADistAtPos(Residue* binderRes);

//         /**
//          * @brief Compute the joint probability from multiple amino acid distributions
//          * 
//          * @param aaDists The inner vector is the probability distribution, the outer vector is each distribution
//          * @return vector<mstreal> 
//          */
//         vector<mstreal> getJointProbability(vector<vector<mstreal>> aaDists);

//     private:
//         augmentedStructure augmentedTarget;
//         augmentedStructure* augmentedBinder = nullptr;

//         mstreal designabilityScore = 0.0;
//         mstreal sequenceCompatibilityScore = 0.0;

//         int numNonDesignable = 0;

//         map<res_t,frameTable*> frameTables;

//         vector<int> interfaceCounts;
//         vector<mstreal> interfaceDesignabilityScores;
//         vector<mstreal> interfaceSequenceCompatibilityScores;
//         vector<int> interfaceTotalNumberOfMatches;
//         vector<mstreal> searchTimes;

//         fstream* binder_info_out = nullptr;
//         fstream* contact_info_out = nullptr;
//         fstream* clash_info_out = nullptr;
//         // fstream* match_info_out = nullptr;
// };

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