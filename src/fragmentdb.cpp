#include "fragmentdb.h"

/* --- --- --- --- --- clusterDBSegments --- --- --- --- --- */

void clusterDBSegments::loadFASSTDB(string fasstDBPath) {
    cout << "Reading fasstDB..." << endl;
    F.readDatabase(fasstDBPath,1);
    cout << "Done reading database..." << endl;
};

void clusterDBSegments::cluster(mstreal RMSDCutoff, mstreal coverage, int nSubsampled) {
    // Reset object
    allSegments.clear();
    clusters.clear();
    if (!clusterRepresentatives.empty()) {
        for (Structure* S: clusterRepresentatives) delete S;
        clusterRepresentatives.clear();
    }

    // Get all segments from the DB
    for (int i = 0; i < F.numTargets(); i++) {
        Structure* target = F.getTarget(i);
        for (int j = 0; j < target->chainSize(); j++) {
            Chain* targetC = &target->getChain(j);
            vector<Atom*> chainAtoms = targetC->getAtoms();
            if (4*targetC->residueSize() != chainAtoms.size()) MstUtils::error("Chain does not have the proper number of atoms. Has "+MstUtils::toString(chainAtoms.size())+", but should have "+MstUtils::toString(targetC->residueSize()));
            for (int k = 0; k < targetC->residueSize() - segmentLength + 1; k++) {
                vector<Atom*> segment(chainAtoms.begin()+(4*(k)),chainAtoms.begin()+(4*(k+segmentLength)));
                allSegments.push_back(segment);
            }
        }
    }
    cout << "Extracted " << allSegments.size() << " segments total from the DB" << endl;

    // Do greedy clustering
    Clusterer Cl;
    cout << "Greedily clustering the segments." << endl;
    cout << "RMSD cutoff: " << RMSDCutoff << endl;
    cout << "coverage: " << coverage << endl;
    cout << "nSubsampled: " << nSubsampled << endl;
    clusters = Cl.greedyCluster(allSegments,RMSDCutoff,nSubsampled,coverage,-1,true);
    cout << "Defined " << clusters.size() << " clusters" << endl;

    // Collect the representatives
    for (int i = 0; i < clusters.size(); i++) {
        Structure* clusterRep = new Structure(allSegments[clusters[i].front()]);
        clusterRep->setName("cluster_"+MstUtils::toString(i)+"_representative");
        CartesianPoint centroid = (AtomPointerVector(clusterRep->getAtoms())).getGeometricCenter();
        for (Atom* A : clusterRep->getAtoms()) A->setCoor(A->getCoor()-centroid);

        clusterRepresentatives.push_back(clusterRep);
    }
}

void clusterDBSegments::writeClusterRepresentativesToBinFile(string path) {
    fstream bin_out;
    MstUtils::openFile(bin_out,path,fstream::out | fstream::binary,"clusterDBSegments::writeClusterRepresentativesToBinFile");
    for (Structure* S: clusterRepresentatives) {
        MstUtils::writeBin(bin_out,'S'); // Structure section
        MstUtils::writeBin(bin_out,S->getName());
        S->writeData(bin_out);
    }
    bin_out.close();
}

void clusterDBSegments::writeClusterRepresentativesToPDBFile(string path, int topN) {
    fstream out;
    MstUtils::openFile(out,path,fstream::out,"clusterDBSegments::writeClusterRepresentativesToPDBFile");
    topN = min(topN,int(clusterRepresentatives.size()));
    for (int i = 0; i < topN; i++) {
        Structure* S = clusterRepresentatives[i];
        out << "HEADER " << S->getName() << endl;
        S->writePDB(out);
    }
    out.close();
}

void clusterDBSegments::writeClusterMembersToPDBFile(string path, int topNclusters, int members) {
    fstream out;
    MstUtils::openFile(out,path,fstream::out,"clusterDBSegments::writeClusterRepresentativesToPDBFile");
    topNclusters = min(topNclusters,int(clusterRepresentatives.size()));
    RMSDCalculator calc;
    for (int i = 0; i < topNclusters; i++) {
        int clusterMembers = min(members,int(clusters[i].size()));
        for (int j = 0; j < clusterMembers; j++) {
            vector<Atom*>& segment = allSegments[clusters[i][j]];
            Structure S(segment);
            calc.align(segment,clusterRepresentatives[i]->getAtoms(),S);
            out << "HEADER " << i << "_" << j << endl;
            S.writePDB(out);
        }
    }
    out.close(); 
}

void clusterDBSegments::writeClusterInfo(string path) {
    fstream out;
    MstUtils::openFile(out,path,fstream::out,"clusterDBSegments::writeClusterInfo");
    out << "cluster,size" << endl;
    for (int i = 0; i < clusters.size(); i++) {
        out << i << "," << clusters[i].size() << endl;
    }
    out.close();
}

vector<Structure*> clusterDBSegments::readClusterRepresentativesFromBinFile(string path) {
    vector<Structure*> segments;
    fstream bin_in;
    MstUtils::openFile(bin_in,path,fstream::in | fstream::binary,"clusterDBSegments::readClusterRepresentativesFromBinFile");
    char sect;
    string name;
    while (bin_in.peek() != EOF) {
        MstUtils::readBin(bin_in, sect);
        if (sect == 'S') {
            Structure* S = new Structure();
            MstUtils::readBin(bin_in,name);
            S->readData(bin_in);
            S->setName(name);
            segments.push_back(S);
        } else {
            MstUtils::error("woops, this should be a structure section","clusterDBSegments::readClusterRepresentativesFromBinFile");
        }
    }
    bin_in.close();
    return segments;
}

/* --- --- --- --- --- segmentGraph --- --- --- --- --- */

void segmentGraph::addNode(Structure* segment) {
    if (allNodesSet.count(segment) != 0) MstUtils::error("Structure has already been added to the graph","segmentGraph::addNode");
    if (name2Structure.count(segment->getName()) != 0) MstUtils::error("Structure with name ("+segment->getName()+") has already been added to the graph","segmentGraph::addNode");
    allNodes.push_back(segment); //graph takes ownership of the Structure
    allNodesSet.insert(segment);
    name2Structure[segment->getName()] = segment;
    adjacencies[segment] = vector<Structure*>();
    revAdjacencies[segment] = vector<Structure*>();
}

void segmentGraph::addNodeDesignability(string nodeName, int designability) {
    if (name2Structure.count(nodeName) == 0) MstUtils::error("Node name: "+MstUtils::toString(nodeName)+" not recognized","segmentGraph::addNodeDesignability");
    name2designability[nodeName] = designability;
}

void segmentGraph::addEdge(Structure* source, Structure* destination) {
    if (allNodesSet.count(source) == 0) MstUtils::error("Source node not found in graph","segmentGraph::addEdge");
    if (allNodesSet.count(destination) == 0) MstUtils::error("Destination node not found in graph","segmentGraph::addEdge");
    // (source segment C-terminus -> destination segment N-terminus)
    adjacencies[source].push_back(destination);
    revAdjacencies[destination].push_back(source);
}

bool segmentGraph::isStructureInGraph(string segmentName) {
    return (name2Structure.count(segmentName) == 1);
}

vector<Structure*> segmentGraph::getNtermNeighbors(string nodeName) {
    if (!isStructureInGraph(nodeName)) MstUtils::error("Node not found in graph","segmentGraph::getNtermNeighbors");
    return adjacencies[name2Structure.at(nodeName)];
}

vector<Structure*> segmentGraph::getCtermNeighbors(string nodeName) {
    if (!isStructureInGraph(nodeName)) MstUtils::error("Node not found in graph","segmentGraph::getCtermNeighbors");
    return revAdjacencies[name2Structure.at(nodeName)];
}

Structure* segmentGraph::getNodeByName(string nodeName) {
    if (!isStructureInGraph(nodeName)) MstUtils::error("Node not found in graph","segmentGraph::getCtermNeighbors");
    return name2Structure[nodeName];
}

int segmentGraph::getNodeDesignability(string nodeName) {
    if (!isStructureInGraph(nodeName)) MstUtils::error("Node not found in graph","segmentGraph::getCtermNeighbors");
    if (name2designability.count(nodeName) == 0) MstUtils::error("Node with name: "+MstUtils::toString(nodeName)+" does not have designability info","segmentGraph::getNodeDesignability");
    return name2designability[nodeName];
}

void segmentGraph::constructGraphFromSegments(vector<Structure*> segments, mstreal RMSDcutoff, vector<int> segmentDesignability) {
    int segmentLength = segments.front()->residueSize();
    if (segmentLength % 2 != 0) MstUtils::error("Segments must have an even number of residues","segmentGraph::constructGraphFromSegments");
    for (Structure* segment : segments) {
        if (segment->residueSize() != segmentLength) MstUtils::error("All segments must have the same number of residues","segmentGraph::constructGraphFromSegments");
        addNode(segment);
    }
    RMSDCalculator calc;
    cout << segments.size() << " segments to compare to find overlaps..." << endl;
    if (segmentDesignability.empty()) segmentDesignability = vector<int>(segments.size(),0);
    if (segmentDesignability.size() != segments.size()) MstUtils::error("Wrong number of segment designability values ("+MstUtils::toString(segmentDesignability.size())+"). Should be "+MstUtils::toString("segmentDesignability"),"segmentGraph::constructGraphFromSegments");
    int idx = 0;
    for (Structure* segment_i : segments) {
        // Find overlaps
        vector<Atom*> segment_i_bb = MiscTools::getBackboneAtoms(*segment_i);
        vector<Atom*> segment_i_nterm(segment_i_bb.begin(),segment_i_bb.begin()+(4*segmentLength/2));
        for (Structure* segment_j : segments) {
            vector<Atom*> segment_j_bb = MiscTools::getBackboneAtoms(*segment_j);
            vector<Atom*> segment_j_cterm(segment_j_bb.begin()+(4*segmentLength/2),segment_j_bb.end());
            // Compare N-terminal half of i to C-terminal half of j
            bool success = true;
            mstreal lRMSD = calc.bestRMSD(segment_i_nterm,segment_j_cterm,false,&success);
            if (!success) MstUtils::error("Unable to calculate RMSD","segmentGraph::constructGraphFromSegments");
            if (lRMSD <= RMSDcutoff) addEdge(segment_i,segment_j);
        }
        // Add designability info
        addNodeDesignability(segment_i->getName(),segmentDesignability[idx]);
        if (idx % 1000 == 0) cout << idx << "..." << endl;
        idx++;
    }

    cout << "Done comparing segment termini" << endl;
}

void segmentGraph::writeGraphFile(string pathPrefix) {
    fstream bin_out;
    MstUtils::openFile(bin_out,pathPrefix+"_segmentGraph.bin",fstream::out | fstream::binary,"segmentGraph::writeGraphFile");
    fstream info_out;
    MstUtils::openFile(info_out,pathPrefix+"_segmentGraph_adjlist.txt",fstream::out,"segmentGraph::writeGraphFile");

    // First write all nodes (structures)
    for (Structure* S: allNodes) {
        MstUtils::writeBin(bin_out,'S'); // Structure section
        MstUtils::writeBin(bin_out,S->getName());
        S->writeData(bin_out);
    }
    // Next write node properties
    for (Structure* S: allNodes) {
        MstUtils::writeBin(bin_out,'D'); // Designability section
        string name = S->getName();
        MstUtils::writeBin(bin_out,name);
        MstUtils::writeBin(bin_out,name2designability[name]);
    }
    // Finally write all edges (structure overlap)
    for (Structure* S: allNodes) {
        MstUtils::writeBin(bin_out,'O'); // Overlap section
        MstUtils::writeBin(bin_out,S->getName());
        info_out << S->getName();
        vector<Structure*> nodeAdjacencies = adjacencies[S];
        MstUtils::writeBin(bin_out,int(nodeAdjacencies.size()));
        for (Structure* adjS : nodeAdjacencies) {
            MstUtils::writeBin(bin_out,adjS->getName());
            info_out << " " << adjS->getName();
        }
        info_out << endl;
    }
    bin_out.close();
    info_out.close();
}

void segmentGraph::readGraphFile(string path) {
    fstream bin_in;
    MstUtils::openFile(bin_in,path,fstream::in | fstream::binary,"segmentGraph::readGraphFile");
    char sect;
    string name;
    int nAdj;
    int count;
    while (bin_in.peek() != EOF) {
        MstUtils::readBin(bin_in, sect);
        if (sect == 'S') {
            Structure* S = new Structure();
            MstUtils::readBin(bin_in,name);
            S->readData(bin_in);
            S->setName(name);
            addNode(S);
        } else if (sect == 'D') {
            MstUtils::readBin(bin_in,name);
            MstUtils::readBin(bin_in,count);
            addNodeDesignability(name,count);
        } else if (sect == 'O') {
            MstUtils::readBin(bin_in,name);
            Structure* sourceStructure = name2Structure.at(name);
            MstUtils::readBin(bin_in,nAdj);
            for (int i = 0; i < nAdj; i++) {
                MstUtils::readBin(bin_in,name);
                Structure* destinationStructure = name2Structure.at(name);
                addEdge(sourceStructure,destinationStructure);
            }
        } else {
            MstUtils::error("Section not recognized","segmentGraph::readGraphFile");
        }
    }
    bin_in.close();

    // If the file doesn't contain designability info for the segments, just set the value to 0 for all
    if (name2designability.empty()) {
        for (Structure* segment : allNodes) addNodeDesignability(segment->getName(),0);
    }
}

/* --- --- --- --- --- sampleSegmentOverlaps --- --- --- --- --- */

vector<pair<Structure*,int>> sampleSegmentOverlaps::getExtensionSegments(Structure* query, extensionDirection terminus) {
    if (!extensionSegments.empty()) {
        for (Structure* segment : extensionSegments) delete segment;
        extensionSegments.clear();
    }
    vector<Structure*> overlappingSegmentsUntransformed;
    bool inGraph = overlapGraph.isStructureInGraph(query->getName());
    // Get extension segments with overlap to the current segment
    if (inGraph) {
        // Segment is in graph, so overlap information is known
        if (terminus == extensionDirection::NTERM) {
            overlappingSegmentsUntransformed = overlapGraph.getNtermNeighbors(query->getName());
        } else if (terminus == extensionDirection::CTERM) {
            overlappingSegmentsUntransformed = overlapGraph.getCtermNeighbors(query->getName());
        } else {
            MstUtils::error("Must select a single terminus","sampleSegmentOverlaps::getExtensionSegments");
        }
    } else {
        // No pre-computed overlap information, will need to search against all segments
        vector<Structure*> segments = overlapGraph.getAllNodes();
        mstreal RMSD;
        for (Structure* segment : segments) {
            RMSD = calculateLeastRMSD(query,segment,terminus,false);
            if (RMSD <= overlapRMSDCutoff) overlappingSegmentsUntransformed.push_back(segment);
        }
    }
    // cout << "Found " << overlappingSegmentsUntransformed.size() << " segments with overlap to the current terminus" << endl;

    // Make a copy of each extension segment and apply transformation to optimally superimpose based on the overlap
    vector<pair<Structure*,int>> result;
    for (Structure* segment : overlappingSegmentsUntransformed) {
        Structure* segmentCopy = new Structure(*segment);
        calculateLeastRMSD(query,segmentCopy,terminus,true);
        extensionSegments.push_back(segmentCopy);
        int numMatches = overlapGraph.getNodeDesignability(segmentCopy->getName());
        result.push_back(pair<Structure*,int>(segmentCopy,numMatches));
    }
    return result;
}

mstreal sampleSegmentOverlaps::calculateLeastRMSD(Structure* query, Structure* candidate, extensionDirection terminus, bool transform) {
    bool succ = false;
    vector<Atom*> queryBBAtom = MiscTools::getBackboneAtoms(*query);
    vector<Atom*> candidateBBAtom = MiscTools::getBackboneAtoms(*candidate);
    vector<Atom*> queryTerminus,candidateTerminus;
    int overlapLength = overlapGraph.getOverlapLength();
    if (overlapLength > query->residueSize()) MstUtils::error("Fragment is too short to be extended","sampleSegmentOverlaps::calculateLeastRMSD");
    // The following lines extract the proper subsets of backbone atoms
    if (terminus == extensionDirection::NTERM) {
        // Extend the N terminus of the query
        queryTerminus = vector<Atom*>(queryBBAtom.begin(),queryBBAtom.begin()+(4*overlapLength));
        candidateTerminus = vector<Atom*>(candidateBBAtom.end()-(4*overlapLength),candidateBBAtom.end());
    } else if (terminus == extensionDirection::CTERM) {
        // Extend the C terminus of the query
        queryTerminus = vector<Atom*>(queryBBAtom.end()-(4*overlapLength),queryBBAtom.end());
        candidateTerminus = vector<Atom*>(candidateBBAtom.begin(),candidateBBAtom.begin()+(4*overlapLength));
    } else {
        MstUtils::error("Must select a single terminus","sampleSegmentOverlaps::calculateLeastRMSD");
    }
    mstreal RMSD = calc.bestRMSD(candidateTerminus,queryTerminus,transform,&succ);
    if (!succ) MstUtils::error("Error calculating overlap RMSD","sampleSegmentOverlaps::calculateLeastRMSD");
    if (transform) calc.applyLastTransformation(candidateBBAtom);
    return RMSD;
}

/* --- --- --- --- --- sampleSegmentOverlaps --- --- --- --- --- */

mstreal searchSegments::findLowestRMSDSegment(Structure* query) {
    Structure queryBBAtoms(MiscTools::getBackboneAtoms(*query));
    if (queryBBAtoms.residueSize() != graph.getSegmentLength()) MstUtils::error("Query does not have the same number of residues as segments in the graph","searchSegments::findLowestRMSDSegment");

    vector<Structure*> segments = graph.getAllNodes();
    mstreal minRMSD = 100.0;
    for (Structure* segment : segments) {
        mstreal newRMSD = calc.bestRMSD(queryBBAtoms.getAtoms(),segment->getAtoms());
        if (newRMSD <= minRMSD) {
            minRMSD = newRMSD;
            matchingSegmentName = segment->getName();
            RMSDToMatch = minRMSD;
        }
    }
    return RMSDToMatch;
}

Structure searchSegments::getClosestMatch() {
    if (matchingSegmentName == "") MstUtils::error("Must find a matching segment first","searchSegments::getClosestMatch");
    return *graph.getNodeByName(matchingSegmentName);
}
int searchSegments::getNumMatchesInDB() {
    if (matchingSegmentName == "") MstUtils::error("Must find a matching segment first","searchSegments::getNumMatchesInDB");
    return graph.getNodeDesignability(matchingSegmentName);
}


