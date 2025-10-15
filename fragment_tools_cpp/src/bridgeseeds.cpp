#include "bridgeseeds.h"

/* --- --- --- --- --- seedBridge --- --- --- --- --- */

seedBridge::seedBridge(int _targetID, Residue* ntermR, Residue* ctermR) {
    targetID = _targetID;
    ntermResIdx = ntermR->getResidueIndex();
    residueLength = ctermR->getResidueIndex() - ntermResIdx - 1;
    CartesianPoint distances = findCADistances(ntermR,ctermR);
    R1CaDistance = distances.getX();
    R2CaDistance = distances.getY();
    R3CaDistance = distances.getZ();
    if ((R1CaDistance == 0)||(R2CaDistance == 0)||(R3CaDistance == 0)) {
        cout << "target: " << targetID << endl;
        cout << "ntermR: " << ntermResIdx << endl;
        cout << "residueL: " << residueLength << endl;
        MstUtils::error("Terminal residue distances are 0.0. This means either the object was not yet initialized, or there was an error when defining the termini");
    }
}

int seedBridge::getUniqueID() {
    if (uniqueID == -1) MstUtils::error("The unique ID has not yet been set");
    return uniqueID;
}

void seedBridge::writeToBin(ostream& ofs) {
    // write structure DB info
    MstUtils::writeBin(ofs,targetID);
    MstUtils::writeBin(ofs,ntermResIdx);
    MstUtils::writeBin(ofs,residueLength);

    // write distance info
    MstUtils::writeBin(ofs,R1CaDistance);
    MstUtils::writeBin(ofs,R2CaDistance);
    MstUtils::writeBin(ofs,R3CaDistance);
}

void seedBridge::readFromBin(istream& ifs) {
    // get structure DB info
    MstUtils::readBin(ifs,targetID);
    MstUtils::readBin(ifs,ntermResIdx);
    MstUtils::readBin(ifs,residueLength);

    // get distance info
    MstUtils::readBin(ifs,R1CaDistance);
    MstUtils::readBin(ifs,R2CaDistance);
    MstUtils::readBin(ifs,R3CaDistance);
}

CartesianPoint seedBridge::getTerminalResidueDistances() {
    if ((R1CaDistance == 0)||(R2CaDistance == 0)||(R3CaDistance == 0)) {
        cout << "target: " << targetID << endl;
        cout << "ntermR: " << ntermResIdx << endl;
        cout << "residueL: " << residueLength << endl;
        MstUtils::error("Terminal residue distances are 0.0. This means either the object was not yet initialized, or there was an error when defining the termini");
    } 
    return CartesianPoint(R1CaDistance,R2CaDistance,R3CaDistance);
 }

CartesianPoint seedBridge::findCADistances(Residue* ntermR, Residue* ctermR) {
    // Select the Ca atoms at the N-terminal end of the bridge
    Residue* ntermR1 = ntermR->iPlusDelta(-2);
    Residue* ntermR2 = ntermR->iPlusDelta(-1);
    Residue* ntermR3 = ntermR;

    if ((ntermR1 == NULL)||(ntermR2 == NULL)||(ntermR3 == NULL)) MstUtils::error("Missing residues at the N-terminus of the bridge","seedBridge::findCADistances");
    CartesianPoint ntermR1_CA = ntermR1->findAtom("CA");
    CartesianPoint ntermR2_CA = ntermR2->findAtom("CA");
    CartesianPoint ntermR3_CA = ntermR3->findAtom("CA");

    // Select the residues at the C-terminal end of the bridge
    Residue* ctermR1 = ctermR;
    Residue* ctermR2 = ctermR->iPlusDelta(1);
    Residue* ctermR3 = ctermR->iPlusDelta(2);
    if ((ctermR1 == NULL)||(ctermR2 == NULL)||(ctermR3 == NULL)) MstUtils::error("Missing residues at the C-terminus of the bridge","seedBridge::findCADistances");
    CartesianPoint ctermR1_CA = ctermR1->findAtom("CA");
    CartesianPoint ctermR2_CA = ctermR2->findAtom("CA");
    CartesianPoint ctermR3_CA = ctermR3->findAtom("CA");

    return CartesianPoint(ntermR1_CA.distance(ctermR1_CA),ntermR2_CA.distance(ctermR2_CA),ntermR3_CA.distance(ctermR3_CA));
}

vector<Atom*> seedBridge::getBridgeTerminusFromStructure(Residue* ntermR, Residue* ctermR, int terminusLength, bool terminusOnly) {
    // if (ntermR->getParent() != ctermR->getParent()) MstUtils:error("N and C-terminal residues must be from the same chain","seedBridge::getBridgeTerminusFromStructure");
    vector<Atom*> terminalResidueAtoms;
    // N-terminal residues
    Residue* R = nullptr;
    for (int i = terminusLength - 1; i >= 0; i--) {
        R = ntermR->iPlusDelta(-i);
        vector<Atom*> bbAtoms = RotamerLibrary::getBackbone(R);
        terminalResidueAtoms.insert(terminalResidueAtoms.end(),bbAtoms.begin(),bbAtoms.end());
    }
    // Central bridge residues
    if (!terminusOnly) {
        R = R->nextResidue();
        while (R != ctermR) {
            vector<Atom*> bbAtoms = RotamerLibrary::getBackbone(R);
            terminalResidueAtoms.insert(terminalResidueAtoms.end(),bbAtoms.begin(),bbAtoms.end());
            R = R->nextResidue();
        }
    }
    // C-terminal residues
    for (int i = 0; i < terminusLength; i++) {
        R = ctermR->iPlusDelta(i);
        vector<Atom*> bbAtoms = RotamerLibrary::getBackbone(R);
        terminalResidueAtoms.insert(terminalResidueAtoms.end(),bbAtoms.begin(),bbAtoms.end());
    }
    return terminalResidueAtoms;
}

/* --- --- --- --- --- seedBridgeDB --- --- --- --- --- */

void seedBridgeDB::buildDBfromStructures(int _maxBridgeLen) {
    maxBridgeLength = _maxBridgeLen;
    allBridges.clear();
    if (structureDB->numTargets() <= 0) MstUtils::error("Must load structures before defining bridge elements","seedBridgeDB::buildDBfromStructures");
    for (int targetID = 0; targetID < structureDB->numTargets(); targetID++) {
        const augmentedStructure* target = &(structureDB->getTarget(targetID));
        for (int chainIDX = 0; chainIDX < target->chainSize(); chainIDX++) {
            Chain* C = &target->getChain(chainIDX);
            for (int bridgeLength = 0; bridgeLength <= maxBridgeLength; bridgeLength++) {
                for (int resID = terminusLength - 1; resID < C->residueSize() - terminusLength - bridgeLength; resID++) {
                    // Now, extract a bridge from residues in chain
                    Residue* ntermR = &C->getResidue(resID);
                    Residue* ctermR = &C->getResidue(resID+bridgeLength+1);
                    shared_ptr<seedBridge> newBridge = make_shared<seedBridge>(targetID,ntermR,ctermR);
                    allBridges.push_back(newBridge);
                }
            }
        }
    }
    cout << "Found " << allBridges.size() << " seed bridges" << endl;
}

void seedBridgeDB::writeDBtoFile(string pathToDBFile) {
    fstream ofs;
    MstUtils::openFile(ofs,pathToDBFile,fstream::out|fstream::binary);

    MstUtils::writeBin(ofs,'V');
    MstUtils::writeBin(ofs, (int)1); // format version
    MstUtils::writeBin(ofs,'L');
    MstUtils::writeBin(ofs, (int)maxBridgeLength);
    for (const shared_ptr<seedBridge>& sB : allBridges) {
        MstUtils::writeBin(ofs,'B');
        sB->writeToBin(ofs);
    }
    MstUtils::writeBin(ofs,'E');
    ofs.close();
}

void seedBridgeDB::readDBfromFile(string pathToDBFile) {
    cout << "Reading bridge binary database from file... filtering out bridges with > " << maxBridgeLength << " residues" << endl;
    fstream ifs;
    MstUtils::openFile(ifs,pathToDBFile,fstream::in|fstream::binary);

    char sect;
    int val;

    MstUtils::readBin(ifs,sect);
    if (sect != 'V') MstUtils::error("File missing version number","seedBridgeDB::readDBfromFile");
    MstUtils::readBin(ifs,val);
    cout << "Loading seed bridge database version " << val << endl;
    MstUtils::readBin(ifs,sect);
    if (sect != 'L') MstUtils::error("File missing bridge length","seedBridgeDB::readDBfromFile");
    MstUtils::readBin(ifs,val);
    cout << "File contains bridges up to length: " << val << endl;
    if (maxBridgeLength < 0) maxBridgeLength = val; // set max to db value if none is already specified
    MstUtils::readBin(ifs,sect);
    int count = 0;
    while (sect == 'B') {
        shared_ptr<seedBridge> sB = make_shared<seedBridge>();
        sB->readFromBin(ifs);
        if ((maxBridgeLength < 0)||(sB->getBridgeLength() <= maxBridgeLength)) {
            // Debugging feature: only load a fraction of protein structures
            // in this case we only keep a seed bridge if it comes from a protein structure that we chose to load
            if ((structureDB!=nullptr)&&(sB->getTargetID() > structureDB->numTargets()-1)) continue;

            sB->setUniqueID(allBridges.size());
            allBridges.push_back(sB);
        }
        // no need to delete, now that we're using shared pointers
        // else {
        //     delete sB;
        // }
        MstUtils::readBin(ifs,sect);
        count++;
        // if (allBridges.size() % 1000 == 0) cout << allBridges.size() << " seed bridges imported so far" << endl;
    }
    ifs.close();
    cout << "Done reading database. Loaded " << allBridges.size() << " bridges" << endl;
}

shared_ptr<seedBridge> seedBridgeDB::getBridge(int uniqueID) {
    if ((uniqueID < 0)||(uniqueID >= allBridges.size())) MstUtils::error("Provied uniqueID falls outside of the bounds of existing seedBridges","seedBridgeDB::getBridge");
    return allBridges[uniqueID];
}

vector<Atom*> seedBridgeDB::getBridgeTerminusFromDB(seedBridge* sB) {
    if (!structureDB) MstUtils::error("Must load protein structure DB in order to extract bridge terminus from DB");
    const Structure& target = structureDB->getTarget(sB->getTargetID());
    return sB->getBridgeTerminusFromStructure(&target.getResidue(sB->getNtermResIdx()),&target.getResidue(sB->getCtermResIdx()),terminusLength,true);
}

Structure seedBridgeDB::getBridgeAndTerminusFromDB(seedBridge* sB) {
    if (!structureDB) MstUtils::error("Must load protein structure DB in order to extract bridge terminus from DB");
    const Structure& target = structureDB->getTarget(sB->getTargetID());
    Structure bridgeAndTerminus = sB->getBridgeTerminusFromStructure(&target.getResidue(sB->getNtermResIdx()),&target.getResidue(sB->getCtermResIdx()),terminusLength,false);
    string newName = MstUtils::toString(sB->getTargetID())+"-"+MstUtils::toString(sB->getNtermResIdx())+"-"+MstUtils::toString(bridgeAndTerminus.residueSize());
    bridgeAndTerminus.setName(newName);
    return bridgeAndTerminus;
}

/* --- --- --- --- --- findSeedBridge --- --- --- --- --- */

void findSeedBridge::reportAPVBoundaries() {
    cout << "xlo: " << PS.getXLow() << endl;
    cout << "xhi: " << PS.getXHigh() << endl;
    cout << "ylo: " << PS.getYLow() << endl;
    cout << "yhi: " << PS.getYHigh() << endl;
    cout << "zlo: " << PS.getZLow() << endl;
    cout << "zhi: " << PS.getZHigh() << endl;
}

void findSeedBridge::setSearchQuery(const shared_ptr<Structure>& nTermSeed, const shared_ptr<Structure>& cTermSeed, int S1CterminusResOffset, int S2NterminusResOffset) {
    // define the query based on distances between the terminal residues
    if ((nTermSeed->residueSize() - S1CterminusResOffset < bridgeData.getTerminusLength())||(cTermSeed->residueSize() - S2NterminusResOffset < bridgeData.getTerminusLength())) 
        MstUtils::error("Provided nTermStructure with n_residues: "+MstUtils::toString(nTermSeed->residueSize())+", offset: "+MstUtils::toString(S1CterminusResOffset)
        +" and cTermStructure with n_residues: "+MstUtils::toString(cTermSeed->residueSize())+", offset: "+MstUtils::toString(S2NterminusResOffset)
        +" not compatible with terminus length: "+MstUtils::toString(bridgeData.getTerminusLength())+")","findSeedBridge::searchSeeds");
    Residue* ntermR = &nTermSeed->getResidue(nTermSeed->residueSize()-1-S1CterminusResOffset);
    Residue* ctermR = &cTermSeed->getResidue(0+S2NterminusResOffset);
    queryDistances = seedBridge::findCADistances(ntermR,ctermR);

    // extract the terminal residues to make a new structure that is used when verifying by superposition
    terminalResidueBBAtoms = seedBridge::getBridgeTerminusFromStructure(ntermR,ctermR,bridgeData.getTerminusLength(),true);
}

void findSeedBridge::setSearchQuery(Structure* nTermSeed, Structure* cTermSeed, int S1CterminusResOffset, int S2NterminusResOffset) {
    
    // define the query based on distances between the terminal residues
    if ((nTermSeed->residueSize() - S1CterminusResOffset < bridgeData.getTerminusLength())||(cTermSeed->residueSize() - S2NterminusResOffset < bridgeData.getTerminusLength())) 
        MstUtils::error("Provided nTermStructure with n_residues: "+MstUtils::toString(nTermSeed->residueSize())+", offset: "+MstUtils::toString(S1CterminusResOffset)
        +" and cTermStructure with n_residues: "+MstUtils::toString(cTermSeed->residueSize())+", offset: "+MstUtils::toString(S2NterminusResOffset)
        +" not compatible with terminus length: "+MstUtils::toString(bridgeData.getTerminusLength())+")","findSeedBridge::searchSeeds");
    Residue* ntermR = &nTermSeed->getResidue(nTermSeed->residueSize()-1-S1CterminusResOffset);
    Residue* ctermR = &cTermSeed->getResidue(0+S2NterminusResOffset);
    queryDistances = seedBridge::findCADistances(ntermR,ctermR);

    // extract the terminal residues to make a new structure that is used when verifying by superposition
    terminalResidueBBAtoms = seedBridge::getBridgeTerminusFromStructure(ntermR,ctermR,bridgeData.getTerminusLength(),true);
}

int findSeedBridge::searchSeedsByCADistance(mstreal distanceCutoff) {
    cout << "distanceCutoff used for searching: " << distanceCutoff << " angstroms" << endl;

    // reset from previous search
    matches.clear();
    verifiedMatches.clear();
    matchesByLength.clear();

    cout << "R1 CA distance: " << queryDistances.getX() << endl;
    cout << "R2 CA distance: " << queryDistances.getY() << endl;
    cout << "R3 CA distance: " << queryDistances.getZ() << endl;

    // If R2 distance is greater than the max, ignore this pair of seeds
    if (queryDistances.getY() > maxSeedDistance) {
        cout << "Seeds are too far apart (" << maxSeedDistance << "), skipping" << endl;
        return 0;
    }

    // search for termini with matching distances
    vector<int> matchingSeedBridgeUniqueIDS = PS.getPointsWithin(queryDistances,0,distanceCutoff,true);
    for (int uniqueID : matchingSeedBridgeUniqueIDS) {
        shared_ptr<seedBridge> bridgeMatch = bridgeData.getBridge(uniqueID);
        matches.push_back(bridgeMatch);
    }
    cout << "Found " << matches.size() << " matching bridges by comparing CA distances" << endl;
    return matches.size();
}

vector<mstreal> findSeedBridge::getBridgeLengthDist() {
    fstream info_out;
    MstUtils::openFile(info_out,"CaDistanceMatches.csv",fstream::out);
    info_out << "numBridgeResidues,numMatches" << endl;
    vector<int> counts(bridgeData.getMaxBridgeLength(),0);
    vector<mstreal> pDist;
    int totalCount = 0;
    for (const shared_ptr<seedBridge>& sB : matches) {
        counts[sB->getBridgeLength()]++;
        totalCount++;
    }
    for (int len = 0; len <= bridgeData.getMaxBridgeLength(); len++) {
        cout << len << "\t" << counts[len] << endl;
        info_out << len << "," << counts[len] << endl;
        pDist.push_back(mstreal(counts[len])/mstreal(totalCount));
    }
    return pDist;
}

void findSeedBridge::loadSeedBridgeDataIntoAPV() {
    cout << "Loading seed bridge terminal residue distances into APV" << endl;
    const vector<shared_ptr<seedBridge>>& allBridges = bridgeData.getAllBridges();
    boundingBox distanceBB(1.0);
    for (const shared_ptr<seedBridge> sB : allBridges) distanceBB.update(sB->getTerminalResidueDistances());
    int Nbuckets = int(ceil(max(max((distanceBB.getXWidth()), (distanceBB.getYWidth())), (distanceBB.getZWidth()))/0.5)); //0.5 Å is the characteristic distance, eq. copied from ProximitySearch MST class
    PS = ProximitySearch(distanceBB.getXMin(),distanceBB.getYMin(),distanceBB.getZMin(),distanceBB.getXMax(),distanceBB.getYMax(),distanceBB.getZMax(),Nbuckets);
    for (const shared_ptr<seedBridge> sB : allBridges) PS.addPoint(sB->getTerminalResidueDistances(),sB->getUniqueID());
    cout << "Done loading into APV." << endl; 
}

int findSeedBridge::verifyMatchesBySuperposition(mstreal RMSDCutoff) {
    if (matches.empty()) MstUtils::error("Must search for matches before verifying","findSeedBridge::verifyMatchesBySuperposition");
    verifiedMatches.clear();
    matchesByLength.clear();

    int count = 0;
    for (const shared_ptr<seedBridge>& sB : matches) {
        vector<Atom*> bridgeTerminus = bridgeData.getBridgeTerminusFromDB(sB.get());
        mstreal RMSDval = calc.bestRMSD(bridgeTerminus,terminalResidueBBAtoms,true);
        // if (count % 100 == 0) cout << "RMSDval: " << RMSDval << endl;
        if (RMSDval <= RMSDCutoff) {
            sB->setRMSD(RMSDval);
            verifiedMatches.push_back(sB);
            // lengthToBridge[sB->getBridgeLength()].push_back(sB);
        }
        count++;
    }
    // Sort by RMSD
     std::sort(verifiedMatches.begin(),verifiedMatches.end());

    // Sort matches by bridge length
    // for (seedBridge* sB : verifiedMatches) {
    //     matchesByLength[sB->getBridgeLength()].push_back(sB);
    // }

    cout << "Verified " << verifiedMatches.size() << " bridges with RMSD < " << RMSDCutoff << " to the query" << endl;
    return verifiedMatches.size();
}

vector<int> findSeedBridge::getVerifiedBridgeLengthDist(string name) {

    vector<int> counts(bridgeData.getMaxBridgeLength()+1,0);
    vector<mstreal> pDist;
    int totalCount = 0;
    for (const shared_ptr<seedBridge>& sB : verifiedMatches) {
        counts[sB->getBridgeLength()]++;
        totalCount++;
    }
    *info_out << name;
    for (int len = 0; len <= bridgeData.getMaxBridgeLength(); len++) {
        // cout << len << "\t" << counts[len] << endl;
        *info_out << "," << counts[len];
        // pDist.push_back(mstreal(counts[len]));
    }
    *info_out << endl;
    return counts;
}

vector<vector<shared_ptr<Structure>>> findSeedBridge::getVerifiedBridgeStructures() {
    vector<vector<shared_ptr<Structure>>> result(bridgeData.getMaxBridgeLength()+1);
    for (const shared_ptr<seedBridge> sB : verifiedMatches) {
        result[sB->getBridgeLength()].push_back(getAlignedBridgeStructure(sB.get()));
    }
    return result;
}

vector<vector<shared_ptr<Structure>>> findSeedBridge::getVerifiedClusteredBridgeStructures(bool verbose) {
    vector<vector<shared_ptr<Structure>>> result, bridge_structures;
    bridge_structures = findSeedBridge::getVerifiedBridgeStructures();
    result.resize(bridge_structures.size());
    for (int bridge_length = 0; bridge_length < bridge_structures.size(); bridge_length++) {
        vector<shared_ptr<Structure>>& lengthNBridges = bridge_structures[bridge_length];
        if (lengthNBridges.empty()) continue;
        vector<shared_ptr<Structure>> clusterReps = clusterBridgeStructures(lengthNBridges,1000,0.75);
        if (verbose) cout << "After clustering, found " << clusterReps.size() << " clusters for bridges of length " << bridge_length << endl;
        result[bridge_length] = clusterReps;
    }
    return result;
}

vector<shared_ptr<Structure>> findSeedBridge::clusterBridgeStructures(const vector<shared_ptr<Structure>>& lengthNBridges, int Nstructures, mstreal clusterRMSD) {
    vector<vector<Atom*>> bridgeAtoms;
    int count = 1;
    for (const shared_ptr<Structure>& s : lengthNBridges) {
        if (count > Nstructures) break;
        bridgeAtoms.push_back(s->getAtoms());
        count++;
    }
    Clusterer cl;
    vector<vector<int>> clusters = cl.greedyCluster(bridgeAtoms,clusterRMSD);
    vector<shared_ptr<Structure>> clusterReps;
    for (int i = 0; i < clusters.size(); i++) clusterReps.push_back(lengthNBridges[clusters[i][0]]);
    return clusterReps;
}

// vector<Structure> findSeedBridge::getRepresentativeForEachLength() {
//     vector<Structure> result;
//     set<int> stored;
//     for (seedBridge* sB : verifiedMatches) {
//         if (stored.count(sB->getBridgeLength()) == 0) {
//             Structure bridgeAndTerminus = bridgeData.getBridgeAndTerminusFromDB(sB);
//             vector<Atom*> bridgeTerminus = bridgeData.getBridgeTerminusFromDB(sB);
//             calc.align(bridgeTerminus,terminalResidueBBAtoms,bridgeAndTerminus);
//             result.push_back(bridgeAndTerminus);
//             stored.insert(sB->getBridgeLength());
//         }
//     }
//     return result;
// }

shared_ptr<Structure> findSeedBridge::getAlignedBridgeStructure(seedBridge* sB) {
    shared_ptr<Structure> bridgeAndTerminus = make_shared<Structure>(bridgeData.getBridgeAndTerminusFromDB(sB));
    vector<Atom*> bridgeTerminus = bridgeData.getBridgeTerminusFromDB(sB);
    calc.align(bridgeTerminus,terminalResidueBBAtoms,*bridgeAndTerminus);
    return bridgeAndTerminus;
}

void findSeedBridge::writeToMultiPDB(string pathToPDB, string bridgeName, int topN) {
    fstream pdb_out;
    MstUtils::openFile(pdb_out,pathToPDB,fstream::out);
    for (int len = 0; len < matchesByLength.size(); len++) {
        int N = 0;
        for (const shared_ptr<seedBridge>& sB : matchesByLength[len]) {
            N++;
            if ((topN != -1)&&(N >= topN)) continue;
            shared_ptr<Structure> bridgeAndTerminus = getAlignedBridgeStructure(sB.get());
            pdb_out << "HEADER    " << bridgeName << "_" << sB->getBridgeLength() << "_" << sB->getUniqueID() << endl;
            bridgeAndTerminus->writePDB(pdb_out);
        }
    }
}

/* --- --- --- --- --- fuseSeedsAndBridge --- --- --- --- --- */

vector<shared_ptr<Structure>> fuseSeedsAndBridge::fuse() {
    all_fused.clear();

    // For each bridge length, cluster, and fuse the representative to the seeds
    for (int lenN = 0; lenN < bridges.size(); lenN++) {
        vector<shared_ptr<Structure>>& lengthNBridges = bridges[lenN];
        if (lengthNBridges.empty()) continue;
        vector<shared_ptr<Structure>> clusterReps = clusterBridgeStructures(lengthNBridges,1000,0.75);
        cout << "After clustering, found " << clusterReps.size() << " clusters for bridges of length " << lenN << endl;
        int bridgeN = 0;
        for (const shared_ptr<Structure>& selectedBridge : clusterReps) {
            bridgeN++;
            // check if clashes with the target
            if ((clashCheck.isStructureSet())&&(clashCheck.checkForClashesToStructure(selectedBridge->getResidues()))) {
                cout << "Bridge " << bridgeN << " clashes and will be discarded" << endl;
                continue;
            }
            string name = seedA->getName()+"-"+MstUtils::toString(seedA->residueSize()-seedAOffset)+"_stitch_loop"+MstUtils::toString(bridgeN)+"-"+MstUtils::toString(lenN)+"_stitch_"+seedB->getName()+"-"+MstUtils::toString(seedB->residueSize()-seedBOffset);
            cout << "Seed with name: " << name << endl;
            *bridge_out << "HEADER    " << name << endl;
            selectedBridge->writePDB(*bridge_out);

            // Fuse seeds + bridge and store
            fusionTopology topology(seedA->residueSize()-seedAOffset+seedB->residueSize()-seedBOffset+lenN);

            vector<Residue*> seedA_res = getFragmentRes(*seedA,0,seedA->residueSize()-seedAOffset);
            topology.addFragment(seedA_res,getFragResIdx(0,seedA->residueSize()-seedAOffset));

            vector<Residue*> seedB_res = getFragmentRes(*seedB,seedBOffset,seedB->residueSize()-seedBOffset);
            topology.addFragment(seedB_res,getFragResIdx(seedA->residueSize()-seedAOffset+lenN,seedB->residueSize()-seedBOffset));

            // add loop last to preserve seed A/B sequence
            topology.addFragment(*selectedBridge,getFragResIdx(seedA->residueSize()-seedAOffset-overlapLength,selectedBridge->residueSize()));

            shared_ptr<Structure> fusedS = make_shared<Structure>(Fuser::fuse(topology,fuserOut,params));
            fusedS->setName(name);
            fusedS->getChain(0).setID("0");
            all_fused.push_back(fusedS);
            bridgeN++;
        }
    }
    return all_fused;
}

void fuseSeedsAndBridge::writeFusedStructuresToPDB() {
    for (const shared_ptr<Structure>& fused : all_fused) {
        *fused_out << "HEADER    " << fused->getName() << endl;
        fused->writePDB(*fused_out);
    }
}

vector<shared_ptr<Structure>> fuseSeedsAndBridge::clusterBridgeStructures(const vector<shared_ptr<Structure>>& lengthNBridges, int Nstructures, mstreal clusterRMSD) {
    vector<vector<Atom*>> bridgeAtoms;
    int count = 1;
    for (const shared_ptr<Structure>& s : lengthNBridges) {
        if (count > Nstructures) break;
        bridgeAtoms.push_back(s->getAtoms());
        count++;
    }
    vector<vector<int>> clusters = cl.greedyCluster(bridgeAtoms,clusterRMSD);
    vector<shared_ptr<Structure>> clusterReps;
    for (int i = 0; i < clusters.size(); i++) clusterReps.push_back(lengthNBridges[clusters[i][0]]);
    return clusterReps;
}
