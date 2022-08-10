#ifndef _CHAINEXTENSION_H
#define _CHAINEXTENSION_H

#include "mstfasst.h"
#include "mstfuser.h"
#include "msttransforms.h"
#include "msttypes.h"

#include <list>

#include "residueframe.h"
#include "fragmentdb.h"
#include "scoreinterface.h"
#include "utilities.h"

// /**
//  * @brief Manages a topology that grows/shrinks by adding segments with partial overlap to one of the termini
//  * 
//  */
// class extensionTopology {
//     public:
//         // The first segment acts as an anchor for the topology
//         extensionTopology(const Structure& segment) {
//             allSegments.emplace_back(0,segment,true);
//         };

//         Structure& getTerminalSegment(extensionDirection terminus);

//         /**
//          * @brief Adds a segment to either N/C termini, with some overhang
//          * s
//          * @param segment The structure to be appended to the topology (passes ownership)
//          * @param extDir The termini that the segment extends from
//          * @param numResExt The number of residues the segment extends the topology by
//          * @param fixed If true, the segment is fixed in space during fusing and cannot be removed
//          */
//         void addSegment(const Structure& segment, extensionDirection extDir, int numResExt, bool fixed = false);
//         bool removeSegment(extensionDirection terminus);

//         /**
//          * @brief Returns a fusion topology given the current state of the object
//          * 
//          * @return fusionTopology 
//          */
//         fusionTopology buildTopology();
//         void writeTopologyToPDB(string namePrefix);

//         Structure& getFusedStructure();
//         void writeFusedStructureToPDB(string name);
//         void writeFusedStructureToPDB(fstream& out, string name, int state);
//         int getResidueSize();

//     private:
//         struct mappedSegment {
//             mappedSegment(int _ntermPos, const Structure& _segment, bool _fixed) : ntermPos(_ntermPos), numRes(_segment.residueSize()), fixed(_fixed), segment(_segment) {};

//             int ctermPos() {
//                 return ntermPos + numRes - 1;
//             }

//             int ntermPos;
//             int numRes;
//             bool fixed;
//             Structure segment;
//         };
//         list<mappedSegment> allSegments;

//         bool fusedStructureCurrent = false;
//         Structure fused;
//         fusionParams fParams;
//         fusionOutput scores;
// };

// class sampleSegmentExtensionsFASST {
// public:
//     sampleSegmentExtensionsFASST(string fasstDBPath);

//     /**
//      * @brief Search for matches to the query that contain extra N- or C-terminal residues 
//      * 
//      * @param querySeg The segment to be extended
//      * @param numResExt The number of residues to try to extend by
//      * @param extDir The direction to extend in
//      * @return int The number of matching segments with viable extensions
//      */
//     int searchForMatchesToFragment(Structure& querySeg, int numResExt, extensionDirection extDir);
    
    
//     Structure* sampleExtensionFragment();
//     int numRemainingSegments() {return indexSampler.numRemaining();}

// private:
//     FASST F;
//     mstreal RMSDcutoff = 0.1;
//     fasstSolutionSet allMatches;

//     // The following information can be used to excise a match from a protein in DB
//     struct matchWithExtension {
//         matchWithExtension(const fasstSolution& _sol, int _alignment, int _totalLength) : sol(_sol), alignment(_alignment), totalLength(_totalLength) {}; 
//         const fasstSolution& sol;
//         // alignment and length will differ from the fasstSolution depending on the extension
//         int alignment = -1;
//         int totalLength = -1;
//     };

//     void attemptToExtend(const fasstSolution& sol, int numResExt, extensionDirection extDir);
//     Structure* getSegment(matchWithExtension match);
    
//     vector<matchWithExtension> viableMatches;
//     fisherYatesShuffle indexSampler;
// };

// class binderChainExtensionFASST {
//     public:
//         binderChainExtensionFASST(extensionDirection _terminus, Structure& _anchorSegment, string fasstDBPath) : terminus(_terminus), sampler(fasstDBPath), anchorSegment(_anchorSegment), extTopology(anchorSegment) {};

//         // Must be set before running
//         void setFixedStructuralContext(Structure& context, const binderScorerParams& sParams);
//         void loadFASSTDB(string fasstDBPath);

//         void setOverlapSegmentLength(int _overlapSegmentLength) {overlapSegmentLength = _overlapSegmentLength;}
//         void setExtensionSegmentLength(int _extensionSegmentLength) {extensionSegmentLength = _extensionSegmentLength;}

//         void extendChainRandom(int totalExtensionLength);

//         void extendChainGreedy(int totalExtensionLength, int NOptionsPerStep = 10);

//         void extendChainBeamSearch(int totalExtensionLength, int NOptionsPerStep = 10, int beamWidth = 10);

//     protected:
//         bool readyToExtend();
//         Structure getTerminalSegment(extensionDirection terminus, bool allowShorter = true);
//         // Structure findExtensionFragmentWithFASST(Structure& query);
//         bool extensionSegmentClash(Structure* extensionFragment);

//         void openTrajectoryFile(string pathToFile);
//         void saveStructureToTrajectoryFile(const Structure& S, string name);

//     private:
//         // General parameters
//         extensionDirection terminus;
//         MstTimer timer;
//         int overlapSegmentLength = 4;
//         int extensionSegmentLength = 4;

//         // Structural fragments
//         Structure anchorSegment;
//         Structure currentWithExtension;
//         augmentedStructure binder;

//         // Fragment sampler
//         sampleSegmentExtensionsFASST sampler;

//         // Fuser
//         extensionTopology extTopology;
//         fusionParams fParams;
//         fusionOutput scores;

//         // Scoring
//         augmentedStructure target;
//         Structure targetBB;
//         clashChecker clashCheck;
//         residueFrameBinderScorer* scorer;

//         fstream traj_out;
// };

// struct scoredBinder {
//     scoredBinder(const Structure& anchorSegment, string _name, mstreal _score = 0.0) : topology(anchorSegment), score(_score), name(_name) {};

//     // bool operator <( const scoredBinder& other) const {
//     //     return score < other.score;
//     // };

//     string getName() const {return name;}
//     augmentedStructure getFusedStructure() {return topology.getFusedStructure();}

//     void setScore(mstreal _score) {score = _score;}

//     extensionTopology topology;
//     mstreal score;
//     string name;
// };

// enum extensionMode {RANDOM, GREEDY, BEAM};

// class binderChainExtension {
//     public:
//         binderChainExtension(Structure& context, string overlapGraphPath, const binderScorerParams& sParams);
//         void setBinderAnchor(Structure& anchor);

//         void extendChain(int lengthToExtend, extensionDirection terminus = NTERM, extensionMode mode = RANDOM, int beamWidth = 5);

//         void coverChainWithExtensionSegments(Chain* chainToCover, extensionDirection terminus, int beamWidth = 10);

//     protected:
//         bool extensionSegmentClash(Structure* extensionFragment);
//         bool binderIntrachainClash(scoredBinder* binder);

//         mstreal scoreBinder(scoredBinder& topology, int segmentDesignability);

//         void openTrajectoryFile(string pathToFile);
//         void closeTrajectoryFile() {traj_out.close();}

//         void openScoreFile(string pathToFile);
//         void closeScoreFile() {score_out.close();}

//     private:
//         MstTimer timer;

//         Structure anchorSegment;
//         string anchorName;

//         // Finding new extension segments
//         sampleSegmentOverlaps overlaps;

//         // Scoring
//         augmentedStructure target;
//         Structure targetBB;
//         clashChecker clashCheck;
//         residueFrameBinderScorer* scorer;

//         // Recording
//         fstream traj_out;
//         fstream score_out;
//         int state_count;
// };

#endif