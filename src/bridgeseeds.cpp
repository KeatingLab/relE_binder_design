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
    vector<Atom*> terminalResidueAtoms;
    // N-terminal residues
    Residue* R = nullptr;
    for (int i = terminusLength - 1; i >= 0; i--) {
        R = ntermR->iPlusDelta(-i);
        vector<Atom*> bbAtoms = RotamerLibrary::getBackbone(R);
        terminalResidueAtoms.insert(terminalResidueAtoms.end(),bbAtoms.begin(),bbAtoms.end());
    }
    // Bridge residues
    if (!terminusOnly) {
        while (R != ctermR) {
            R = R->nextResidue();
            vector<Atom*> bbAtoms = RotamerLibrary::getBackbone(R);
            terminalResidueAtoms.insert(terminalResidueAtoms.end(),bbAtoms.begin(),bbAtoms.end());
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
    for (seedBridge* sB : allBridges) delete sB;
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
                    seedBridge* newBridge = new seedBridge(targetID,ntermR,ctermR);
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
    for (seedBridge* sB : allBridges) {
        MstUtils::writeBin(ofs,'B');
        sB->writeToBin(ofs);
    }
    MstUtils::writeBin(ofs,'E');
    ofs.close();
}

void seedBridgeDB::readDBfromFile(string pathToDBFile) {
    cout << "Reading bridge binary database from file..." << endl;
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
    maxBridgeLength = val;
    MstUtils::readBin(ifs,sect);
    while (sect == 'B') {
        seedBridge* sB = new seedBridge;
        sB->readFromBin(ifs);
        sB->setUniqueID(allBridges.size());
        allBridges.push_back(sB);
        MstUtils::readBin(ifs,sect);
        // if (allBridges.size() % 1000 == 0) cout << allBridges.size() << " seed bridges imported so far" << endl;
    }
    ifs.close();
    cout << "Done reading database. Loaded " << allBridges.size() << " bridges" << endl;
}

seedBridge* seedBridgeDB::getBridge(int uniqueID) {
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
    return sB->getBridgeTerminusFromStructure(&target.getResidue(sB->getNtermResIdx()),&target.getResidue(sB->getCtermResIdx()),terminusLength,false);
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

void findSeedBridge::setSearchQuery(Structure* S1, Structure* S2) {
    nTermSeed = S1;
    cTermSeed = S2;
    
    // define the query based on distances between the terminal residues
    if ((S1->residueSize() < bridgeData.getTerminusLength())||(S2->residueSize() < bridgeData.getTerminusLength())) MstUtils::error("At least one of the provided structures are shorter than the terminus length ("+MstUtils::toString(bridgeData.getTerminusLength())+")","findSeedBridge::searchSeeds");
    Residue* ntermR = &S1->getResidue(S1->residueSize()-1);
    Residue* ctermR = &S2->getResidue(0);
    queryDistances = seedBridge::findCADistances(ntermR,ctermR);

    // extract the terminal residues to make a new structure that is used when verifying by superposition
    terminalResidueBBAtoms = seedBridge::getBridgeTerminusFromStructure(ntermR,ctermR,bridgeData.getTerminusLength(),true);
}

int findSeedBridge::searchSeedsByCADistance(mstreal distanceCutoff) {
    cout << "distanceCutoff used for searching: " << distanceCutoff << " Å" << endl;

    // reset from previous search
    matches.clear();
    verifiedMatches.clear();

    cout << "R1 CA distance: " << queryDistances.getX() << endl;
    cout << "R2 CA distance: " << queryDistances.getY() << endl;
    cout << "R3 CA distance: " << queryDistances.getZ() << endl;

    // search for termini with matching distances
    vector<int> matchingSeedBridgeUniqueIDS = PS.getPointsWithin(queryDistances,0,distanceCutoff,true);
    for (int uniqueID : matchingSeedBridgeUniqueIDS) {
        seedBridge* bridgeMatch = bridgeData.getBridge(uniqueID);
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
    for (seedBridge* sB : matches) {
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
    const vector<seedBridge*>& allBridges = bridgeData.getAllBridges();
    boundingBox distanceBB(1.0);
    for (seedBridge* sB : allBridges) distanceBB.update(sB->getTerminalResidueDistances());
    int Nbuckets = int(ceil(max(max((distanceBB.getXWidth()), (distanceBB.getYWidth())), (distanceBB.getZWidth()))/0.5)); //0.5 Å is the characteristic distance, eq. copied from ProximitySearch MST class
    PS = ProximitySearch(distanceBB.getXMin(),distanceBB.getYMin(),distanceBB.getZMin(),distanceBB.getXMax(),distanceBB.getYMax(),distanceBB.getZMax(),Nbuckets);
    for (seedBridge* sB : allBridges) PS.addPoint(sB->getTerminalResidueDistances(),sB->getUniqueID());
    cout << "Done loading into APV." << endl; 
}

int findSeedBridge::verifyMatchesBySuperposition(mstreal RMSDCutoff) {
    if (matches.empty()) MstUtils::error("Must search for matches before verifying","findSeedBridge::verifyMatchesBySuperposition");
    for (seedBridge* sB : matches) {
        vector<Atom*> bridgeTerminus = bridgeData.getBridgeTerminusFromDB(sB);
        mstreal RMSDval = calc.bestRMSD(bridgeTerminus,terminalResidueBBAtoms,true);
        if (RMSDval <= RMSDCutoff) {
            verifiedMatches.push_back(sB);
        }
    }
    cout << "Verified " << verifiedMatches.size() << " bridges with RMSD < " << RMSDCutoff << " to the query" << endl;
    return verifiedMatches.size();
}

vector<mstreal> findSeedBridge::getVerifiedBridgeLengthDist() {
    fstream info_out;
    MstUtils::openFile(info_out,"RMSDMatches.csv",fstream::out);
    info_out << "numBridgeResidues,numMatches" << endl;

    vector<int> counts(bridgeData.getMaxBridgeLength(),0);
    vector<mstreal> pDist;
    int totalCount = 0;
    for (seedBridge* sB : verifiedMatches) {
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

vector<Structure> findSeedBridge::getVerifiedBridgeStructures(int bridgeLength) {
    vector<Structure> result;
    for (seedBridge* sB : verifiedMatches) {
        if ((bridgeLength == -1)||(sB->getBridgeLength() == bridgeLength)) {
            Structure bridgeAndTerminus = bridgeData.getBridgeAndTerminusFromDB(sB);
            vector<Atom*> bridgeTerminus = bridgeData.getBridgeTerminusFromDB(sB);
            calc.align(bridgeTerminus,terminalResidueBBAtoms,bridgeAndTerminus);
            result.push_back(bridgeAndTerminus);
        }
    }
    return result;
}

vector<Structure> findSeedBridge::getRepresentativeForEachLength() {
    vector<Structure> result;
    set<int> stored;
    for (seedBridge* sB : verifiedMatches) {
        if (stored.count(sB->getBridgeLength()) == 0) {
            Structure bridgeAndTerminus = bridgeData.getBridgeAndTerminusFromDB(sB);
            vector<Atom*> bridgeTerminus = bridgeData.getBridgeTerminusFromDB(sB);
            calc.align(bridgeTerminus,terminalResidueBBAtoms,bridgeAndTerminus);
            result.push_back(bridgeAndTerminus);
            stored.insert(sB->getBridgeLength());
        }
    }
    return result;
}