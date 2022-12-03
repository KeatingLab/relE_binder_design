#ifndef _FRAGMENTDB_H
#define _FRAGMENTDB_H

#include "mstfasst.h"
#include "mstsequence.h"
#include "msttypes.h"

#include "utilities.h"

enum extensionDirection {NTERM, CTERM, EITHER};

class clusterDBSegments {
    public:
        clusterDBSegments(int _segmentLength = 8) {
            if (_segmentLength % 2 != 0) MstUtils::error("Segment length must be even","clusterDBSegments::clusterDBSegments");
            segmentLength = _segmentLength;
        }
        ~clusterDBSegments() {
            for (Structure* S: clusterRepresentatives) delete S;
        }

        void loadFASSTDB(string fasstDBPath);

        void cluster(mstreal RMSDCutoff = 0.25, mstreal coverage = 0.9, int nSubsampled = 1000);

        void writeClusterRepresentativesToBinFile(string path);
        void writeClusterRepresentativesToPDBFile(string path, int topN);
        void writeClusterMembersToPDBFile(string path, int topNclusters, int members);
        void writeClusterInfo(string path);

        static vector<Structure*> readClusterRepresentativesFromBinFile(string path);

    private:
        int segmentLength;
        FASST F;

        vector<vector<Atom*>> allSegments;
        vector<vector<int>> clusters;
        vector<Structure*> clusterRepresentatives;
};

class segmentGraph {
    public:
        segmentGraph() {;}
        segmentGraph(string overlapGraphPath) {readGraphFile(overlapGraphPath);}
        ~segmentGraph() {
            for (Structure* node : allNodes) delete node;
        }

        void addNode(Structure* segment);
        void addNodeDesignability(string nodeName, int designability);
        void addEdge(Structure* source, Structure* destination);

        bool isStructureInGraph(string segmentName);

        vector<Structure*> getNtermNeighbors(string nodeName);
        vector<Structure*> getCtermNeighbors(string nodeName);
        vector<Structure*> getAllNodes() {return allNodes;}
        Structure* getNodeByName(string nodeName);
        int getNodeDesignability(string nodeName);
        int getSegmentLength() {return allNodes.front()->residueSize();}
        int getOverlapLength() {return allNodes.front()->residueSize()/2;}

        void constructGraphFromSegments(vector<Structure*> segments, mstreal RMSDcutoff = 0.25, vector<int> segmentDesignability = {});

        void writeGraphFile(string pathPrefix);
        void readGraphFile(string path);

    private:
        vector<Structure*> allNodes;
        unordered_set<Structure*> allNodesSet;
        unordered_map<string,int> name2designability;
        unordered_map<string,Structure*> name2Structure;
        unordered_map<Structure*,vector<Structure*>> adjacencies;
        unordered_map<Structure*,vector<Structure*>> revAdjacencies;
};

class sampleSegmentOverlaps {
    public: 
        sampleSegmentOverlaps(string overlapGraphPath, mstreal _overlapRMSDCutoff = 0.5) : overlapGraph(overlapGraphPath), overlapRMSDCutoff(_overlapRMSDCutoff) {;}

        ~sampleSegmentOverlaps() {
            for (Structure* segment : extensionSegments) delete segment;
        }

        int getOverlapLength() {return overlapGraph.getOverlapLength();}

        /**
         * @brief Given some backbone segment, find other segments that overlap to the selected terminus and could be fused to grow the structure
         * 
         * @param query The backbone segment that is to be extended
         * @param Nterm If true, find segments that overlap to the N-terminus, if false, the opposite.
         * @return vector<pair<Structure*,int>> Pointers to segment structures and the number of matches (memory managed by this class)
         */
        vector<pair<Structure*,int>> getExtensionSegments(Structure* query, extensionDirection terminus = NTERM);

    protected:
        mstreal calculateLeastRMSD(Structure* query, Structure* candidate, extensionDirection terminus = NTERM, bool transform = false);

    private:
        segmentGraph overlapGraph;
        mstreal overlapRMSDCutoff;

        RMSDCalculator calc;

        vector<Structure*> extensionSegments;
};

class searchSegments {
    public:
        searchSegments(string overlapGraphPath) : graph(overlapGraphPath) {;}

        mstreal findLowestRMSDSegment(Structure* query);

        int getSegLen() {return graph.getSegmentLength();}
        Structure getClosestMatch();
        int getNumMatchesInDB();

    private:
        segmentGraph graph;

        string matchingSegmentName = "";
        mstreal RMSDToMatch = -1.0;

        RMSDCalculator calc;
};

#endif