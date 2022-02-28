#include "scoreinterface.h"

/* --- --- --- --- --- interfaceSearch --- --- --- --- --- */

interfaceSearch::interfaceSearch(string frameDBPath, mstreal _posCut, mstreal _oriCut) {
    if (_posCut > 0) posCut = _posCut;
    if (_oriCut > 0) oriCut = _oriCut;
    aaTypes = SeqToolsExtension::getAANames();
    frameDB DB(frameDBPath);

    timer.start();

    vector<mobileFrame*> loaded = DB.loadAllFrames();
    cout << "Loaded " << loaded.size() << " frames in total." << endl;

    for (string aa : aaTypes) {
        cout << "Loading frames with amino acid type: " << aa << endl;
        vector<mobileFrame*> currentFrames;
        res_t current_aa = SeqTools::aaToIdx(aa);
        for (mobileFrame* frame : loaded) {
            if (frame->getResIIndex() == current_aa) {
                currentFrames.push_back(frame);
            }
        }
        boundingBox bbox;
        bool verbose = true;
        bbox.update(currentFrames);
        bbox.printBounds();
        frameTable* fTable = new frameTable(bbox, 0.25, 36);
        for (mobileFrame* frame : currentFrames) fTable->insertFrame(frame);
        frameTables[current_aa] = fTable;
    }

    timer.stop();
    cout << "It took " << timer.getDuration() << " seconds to load all frames into hash tables" << endl;
}

void interfaceSearch::searchInterface(alignInteractingFrames* alignF) {
    // reset all variables
    matchingFrames.clear();
    searchTimes.clear();

    // score each interface residue pair
    for (auto resPair : interfaceResidues) {
        searchResiduePair(resPair,alignF);
    }
}

void interfaceSearch::searchResiduePair(pair<Residue*,Residue*> resPair, alignInteractingFrames* alignF) {
    string resAAName = resPair.first->getName();
    if (aaTypes.find(resAAName) == aaTypes.end()) {
        // MstUtils::error("Amino acid  type:"+resAAName+" not a canonical amino acid");
        cout << "Warning: amino acid type: " << resAAName << " is not a canonical amino acid type. Skipping contact" << endl;
        matchingFrames.push_back(set<mobileFrame*>());
        return;
    }
    res_t aaType = SeqTools::aaToIdx(resAAName);
    
    timer.start();
    residueFrame* query = defineQuery(resPair.first,resPair.second);
    set<mobileFrame*> matches = frameTables[aaType]->findSimilarFrames(query,posCut,oriCut);
    timer.stop();

    mstreal duration = timer.getDuration(MstTimer::timeUnits::msec);
    cout << "Found " << matches.size() << " matches in " << duration << " ms" << endl;
    searchTimes.push_back(duration);

    if (alignF != nullptr) {
        set<mobileFrame*> filteredMatches;
        int count = 0;
        for (mobileFrame* frame : matches) {
            bool isHomolog = alignF->isQueryHomologousToMatchInDB(resPair.first,resPair.second,frame);
            if (!isHomolog) filteredMatches.insert(frame);
            count++;
        }
        cout << "After filtering, we have " << filteredMatches.size() << " matches left" << endl;
        matchingFrames.push_back(filteredMatches);
    } else {
        matchingFrames.push_back(matches);
    }

    delete query;
}

void interfaceSearch::writeContactScoresToFile(string name, bool append) {
    if (contact_info_out == nullptr) {
        contact_info_out = new fstream;
        string filename = name + "_contactScores.tsv";
        MstUtils::openFile(*contact_info_out,filename,fstream::out);
        *contact_info_out << "structure_name\ttargetRes\ttargetResName\tbinderRes\tbinderResName\tnumMatches\tnumMatchesWithNativeAA\tdistance\tsearchTime" << endl;
    }
    for (int i = 0; i < interfaceResidues.size(); i++) {
        Residue* targetR = interfaceResidues[i].first;
        Residue* binderR = interfaceResidues[i].second;
        mstreal distanceCA = (targetR->findAtom("CA")->getCoor() - binderR->findAtom("CA")->getCoor()).norm();
        *contact_info_out << binderR->getStructure()->getName() << "\t";
        *contact_info_out << targetR->getChainID() << targetR->getNum() << "\t" << targetR->getName() << "\t";
        *contact_info_out << binderR->getChainID() << binderR->getNum() << "\t" << binderR->getName()<< "\t";
        *contact_info_out << matchingFrames[i].size() << "\t";
        int nativeAAMatch = 0;
        for (mobileFrame* frame : matchingFrames[i]) if (SeqTools::aaToIdx(binderR->getName()) == frame->getResJIndex()) nativeAAMatch++;
        *contact_info_out << nativeAAMatch << "\t";
        *contact_info_out << distanceCA << "\t";
        *contact_info_out << searchTimes[i];
        *contact_info_out << endl;
    }
}

void interfaceSearch::writeContactMatchesToFile(string name, bool append) {
    if (match_info_out == nullptr) {
        match_info_out = new fstream;
        string filename = name + "_contactMatches.tsv";
        MstUtils::openFile(*match_info_out,filename,fstream::out);
        *match_info_out << "structure_name\ttargetRes\ttargetResName\tbinderRes\tbinderResName\tmatchTarget\tmatchResI\tmatchResJ\tmatchResJAA" << endl;
    }
    for (int i = 0; i < interfaceResidues.size(); i++) {
        Residue* targetR = interfaceResidues[i].first;
        Residue* binderR = interfaceResidues[i].second;
        for (mobileFrame* frame : matchingFrames[i]) {
            *match_info_out << binderR->getStructure()->getName() << "\t";
            *match_info_out << targetR->getChainID() << targetR->getNum() << "\t" << targetR->getName() << "\t";
            *match_info_out << binderR->getChainID() << binderR->getNum() << "\t" << binderR->getName()<< "\t";
            *match_info_out << frame->getTarget() << "\t" << frame->getResI() << "\t" << frame->getResJ() << "\t";
            *match_info_out << SeqTools::idxToTriple(frame->getResJIndex());
            *match_info_out << endl;
        }
    }
}

void interfaceSearch::writeMatchStructures(string name, alignInteractingFrames& alignF) {
    fstream pdb_out;
    MstUtils::openFile(pdb_out, name, fstream::out, "interfaceSearch::writeMatchStructures");

    // transform interacting residues
    for (int i = 0; i < interfaceResidues.size(); i++) {
        residueFrame* refFrame = complex->getResidueFrame(interfaceResidues[i].first->getResidueIndex());
        Transform globalToRef = TransformFactory::switchFrames(*refFrame,Frame());

        int count = 0;
        for (mobileFrame* mFrame : matchingFrames[i]) {
            Structure *interactingResiduePair = alignF.getAlignedInteractingRes(mFrame);
            globalToRef.apply(interactingResiduePair);
            pdb_out << "HEADER    " << interactingResiduePair->getName() << "_match" << count << endl;
            interactingResiduePair->writePDB(pdb_out);
            delete interactingResiduePair;
            count++;
        }
    }
    pdb_out.close();
}

void interfaceSearch::writeContactPropertyToFile(string name) {
    fstream conts_out;
    string filename = name + "_contactProperty.tsv";
    MstUtils::openFile(conts_out, filename, fstream::out, "interfaceSearch::writeContactPropertyToFile");

    conts_out << "RiChain\tRiNum\tRjChain\tRjNum\tVal" << endl;
    for (int i = 0; i < interfaceResidues.size(); i++) {
        Residue* Ri = interfaceResidues[i].first;
        Residue* Rj = interfaceResidues[i].second;
        conts_out << Ri->getChainID() << "\t" << Ri->getNum() << "\t";
        conts_out << Rj->getChainID() << "\t" << Rj->getNum() << "\t";
        mstreal val = (matchingFrames[i].empty()) ? 1 : -1; //negative score is considered more favorable, so -1 means covered
        conts_out << val;
        conts_out << endl;
    }
}

/* --- --- --- --- --- binderScorer --- --- --- --- --- */

binderScorer::binderScorer(const binderScorerParams& params, augmentedStructure& _target) : target(_target), posCut(params.posCut), oriCut(params.oriCut)  {
    aaTypes = SeqToolsExtension::getAANames();
    targetName = MstSys::splitPath(target.getName(),1);
    MiscTools::extractBackboneFromStructure(target,targetBackbone);
    prepareVoxelGrids(params.frameDBPath);
}

binderScorer::binderScorer(const binderScorerParams& params, augmentedStructure& complex, string binderChainIDsString, string targetChainIDsString) : posCut(params.posCut), oriCut(params.oriCut) {
    aaTypes = SeqToolsExtension::getAANames();

    vector<string> targetChainIDs = MstUtils::split(targetChainIDsString,"_");
    cout << "target chain IDs:";
    for (string s : targetChainIDs) cout << "\t" << s << endl;

    vector<string> binderChainIDs = MstUtils::split(binderChainIDsString,"_");
    cout << "binder chain IDs:";
    for (string s : binderChainIDs) cout << "\t" << s << endl;

    // extract target chains from complex
    vector<Residue*> targetResidues;
    for (string targetChainID : targetChainIDs) {
        vector<Residue*> targetChainResidues = complex.getChainByID(targetChainID)->getResidues();
        targetResidues.insert(targetResidues.end(),targetChainResidues.begin(),targetChainResidues.end());
    }
    target = Structure(targetResidues);
    targetName = MstSys::splitPath(complex.getName(),1);

    // extract binder chains from complex
    vector<Residue*> binderResidues;
    for (string binderChainID : binderChainIDs) {
        vector<Residue*> binderChainResidues = complex.getChainByID(binderChainID)->getResidues();
        binderResidues.insert(binderResidues.end(),binderChainResidues.begin(),binderChainResidues.end());
    }
    binder = new augmentedStructure(binderResidues);
    binder->setName(targetName);
    complexMode = true; // will need to remember to delete this chain
    
    MiscTools::extractBackboneFromStructure(target,targetBackbone);
    prepareVoxelGrids(params.frameDBPath);

    defineInterfaceUsingVDWContacts();

    setFrameProbabilityTables();
}

void binderScorer::setTargetBindingSiteResidues(vector<Residue*> sel) {
    targetResidues.clear();
    for (Residue* R : sel) targetResidues.insert(R);
    cout << "Set " << targetResidues.size() << " target residues" << endl;
    setFrameProbabilityTables();
}

void binderScorer::defineTargetBindingSiteResiduesByrSASA(mstreal relSASAthreshold) {
    targetResidues.clear();
    sasaCalculator calc(target);
    bool relative = true;
    map<Residue*,mstreal> relSASA = calc.getResidueSASA(relative);
    for (auto it : relSASA) {
        if (it.second >= relSASAthreshold) targetResidues.insert(it.first);
    }
    cout << "Set " << targetResidues.size() << " target residues with relSASA threshold " << relSASAthreshold << endl;
    setFrameProbabilityTables();
}

void binderScorer::defineInterfaceByPotentialContacts() {
    if (targetResidues.empty()) MstUtils::error("Must define target residues before interface");
    vector<Residue*> targetBSRes(targetResidues.begin(),targetResidues.end());
    potentialContacts pC(targetBSRes,getBinderResidues());
    interfaceResidues = pC.getContacts(true);
    cout << "Defined " << interfaceResidues.size() << " potential contacts" << endl;
    // check that there are no target residues in the contacts that remain undefined
    for (auto pair : interfaceResidues) {
        Residue* R = pair.first;
        if (targetResidues.find(R) == targetResidues.end()) MstUtils::error("Could not find target residue: "+R->getChainID()+MstUtils::toString(R->getNum()));
    }
}

void binderScorer::scoreInterface() {
    if (interfaceResidues.empty()) {
        cout << "No contacts to score" << endl;
        return;
        // MstUtils::error("Must define an interface before scoring","interfaceScorer::scoreInterface");
    }
    // reset all variables
    probabilities.clear();
    countAndNorm.clear();
    searchTimes.clear();

    // score each interface residue pair
    for (auto resPair : interfaceResidues) {
        scoreResiduePair(resPair);
    }
}

void binderScorer::scoreResiduePair(pair<Residue*,Residue*> resPair) {
    string resAAName = resPair.first->getName();
    if (aaTypes.find(resAAName) == aaTypes.end()) {
        // MstUtils::error("Amino acid  type:"+resAAName+" not a canonical amino acid");
        cout << "Warning: amino acid type: " << resAAName << " is not a canonical amino acid type. Skipping contact" << endl;
        probabilities.push_back(-1.0);
        return;
    }
    res_t aaType = SeqTools::aaToIdx(resAAName);
    Residue* targetRes = resPair.first;
    
    timer.start();
    residueFrame* query = defineQuery(resPair.first,resPair.second);
    frameTables[aaType]->setSearchParams(posCut,oriCut);
    pair<int,int> vals = frameTables[aaType]->findInteractingFrameProbability(query,targetRes->getResidueIndex());
    timer.stop();

    mstreal duration = timer.getDuration(MstTimer::timeUnits::msec);
    cout << "Found " << vals.first << " matches out of " << vals.second << " in " << duration << " ms" << endl;
    searchTimes.push_back(duration);

    countAndNorm.push_back(vals);
    probabilities.push_back(logProb(vals.first,vals.second));

    delete query;
}

void binderScorer::writeContactScoresToFile(bool append) {
    if (contact_info_out == nullptr) {
        contact_info_out = new fstream;
        string filename = targetName + "_contactScores.tsv";
        MstUtils::openFile(*contact_info_out,filename,fstream::out);
        *contact_info_out << "targetName\tbinderName\ttargetRes\ttargetResName\tbinderRes\tbinderResName\tprobability\tcounts\tnormalizingfactor\tdistance\tsearchTime_ms" << endl;
    }
    for (int i = 0; i < interfaceResidues.size(); i++) {
        Residue* targetR = interfaceResidues[i].first;
        Residue* binderR = interfaceResidues[i].second;
        mstreal distanceCA = (targetR->findAtom("CA")->getCoor() - binderR->findAtom("CA")->getCoor()).norm();
        *contact_info_out << targetName << "\t";
        *contact_info_out << binderR->getStructure()->getName() << "\t";
        *contact_info_out << targetR->getChainID() << targetR->getNum() << "\t" << targetR->getName() << "\t";
        *contact_info_out << binderR->getChainID() << binderR->getNum() << "\t" << binderR->getName()<< "\t";
        *contact_info_out << probabilities[i] << "\t";
        *contact_info_out << countAndNorm[i].first << "\t";
        *contact_info_out << countAndNorm[i].second << "\t";
        *contact_info_out << distanceCA << "\t";
        *contact_info_out << searchTimes[i];
        *contact_info_out << endl;
    }
}

void binderScorer::writeContactPropertyToFile(string dirPath) {
    fstream conts_out;
    if (dirPath != "") dirPath = dirPath + "/";
    string filename = dirPath + binder->getName() + "_contactScoresSimple.tsv";
    cout << "filename: " << filename << endl;
    MstUtils::openFile(conts_out, filename, fstream::out, "interfaceScorer::writeContactPropertyToFile");

    for (int i = 0; i < interfaceResidues.size(); i++) {
        Residue* Ri = interfaceResidues[i].first;
        Residue* Rj = interfaceResidues[i].second;
        mstreal val = probabilities[i];
        conts_out << Ri->getChainID() << Ri->getNum() << "\t";
        conts_out << Rj->getChainID() << Rj->getNum() << "\t";
        conts_out << val;
        conts_out << endl;
    }
}

void binderScorer::writePSSM(string dirPath) {
    fstream pssm_out;
    if (dirPath != "") dirPath = dirPath + "/";
    string filename = dirPath + binder->getName() + "_binderPSSM.csv";
    cout << "filename: " << filename << endl;
    MstUtils::openFile(pssm_out, filename, fstream::out, "interfaceScorer::writePSSM");
    
    // write header
    set<string> aaNames = SeqToolsExtension::getAANames();
    vector<string> aaNamesVec(aaNames.begin(),aaNames.end());
    pssm_out << "position," << MstUtils::join(",",aaNamesVec) << endl;    

    for (Residue* binderR : binder->getResidues()) {
        // if residue has potential contacts, we can estimate the amino acid probability distribution
        vector<mstreal> aaDist = getAADistAtPos(binderR);
        if (aaDist.empty()) continue;
        pssm_out << binderR->getResidueIndexInChain()+1; // seeds are 1-indexed
        for (mstreal val : aaDist) {
            pssm_out << "," << val;
        }
        pssm_out << endl;
    }
    pssm_out.close();
}

void binderScorer::prepareVoxelGrids(string frameDBPath) {
    timer.start();

    frameDB DB(frameDBPath);
    vector<mobileFrame*> loaded = DB.loadAllFrames();
    cout << "Loaded " << loaded.size() << " frames in total." << endl;

    for (string aa : aaTypes) {
        cout << "Loading frames with amino acid type: " << aa << endl;
        vector<mobileFrame*> currentFrames;
        res_t current_aa = SeqTools::aaToIdx(aa);
        for (mobileFrame* frame : loaded) {
            if (frame->getResIIndex() == current_aa) {
                currentFrames.push_back(frame);
            }
        }
        boundingBox bbox;
        bool verbose = true;
        bbox.update(currentFrames);
        bbox.printBounds();
        frameProbability* fTable = new frameProbability(bbox, posCut/mstreal(2), 360/int(oriCut/mstreal(2)));
        for (mobileFrame* frame : currentFrames) fTable->insertFrame(frame);
        frameTables[current_aa] = fTable;
    }

    timer.stop();
    cout << "It took " << timer.getDuration() << " seconds to load all frames into hash tables" << endl;
}

void binderScorer::defineInterfaceUsingVDWContacts() {
    if (!complexMode) MstUtils::error("Can only define interface using van der Waals contacts when scoring a full-atom complex");
    vdwContacts C(target.getResidues(),binder->getResidues());
    interfaceResidues = C.getInteractingRes();
    for (auto pair : interfaceResidues) targetResidues.insert(pair.first);

    cout << interfaceResidues.size() << " interface contacts when defined by van der Waals interactions" << endl;
}

void binderScorer::setFrameProbabilityTables() {
    // for each target residue with the same amino acid identity, compute the occupied volume
    for (string aa : aaTypes) {
        cout << "Set target residues with identity: " << aa << endl;
        res_t current_aa = SeqTools::aaToIdx(aa);
        frameProbability* fTable = frameTables.at(current_aa);
        for (Residue* targetRes : targetResidues) {
            if (SeqTools::aaToIdx(targetRes->getName()) == current_aa) {
                if (fTable->isTargetResidueDefined(targetRes->getResidueIndex())) continue;
                fTable->setTargetResidue(targetRes->getResidueIndex(),targetBackbone);
            }
        }
    }
}

mstreal binderScorer::logProb(mstreal numerator, mstreal denominator, mstreal pseudocount) {
    return log(numerator+pseudocount) - log(denominator+pseudocount);
}

vector<mstreal> binderScorer::getAADistAtPos(Residue* binderRes) {
    vector<mstreal> aaDist;

    // collect the interacting pairs
    vector<pair<Residue*,Residue*>> contactsWithBinderR;
    for (auto pair : interfaceResidues) {
        if (pair.second == binderRes) {
            contactsWithBinderR.push_back(pair);
        }
    }
    if (contactsWithBinderR.empty()) return aaDist;

    // estimate the 'diagonal' of the joint distribution 
    // e.g. what is the probability of some amino acid at this binder residue, given
    // each of the target residues (assuming independence)
    vector<vector<mstreal>> allContactAADist;
    for (pair<Residue*,Residue*> resPair : contactsWithBinderR) {
        string resAAName = resPair.first->getName();
        if (aaTypes.find(resAAName) == aaTypes.end()) {
            cout << "Warning: amino acid type: " << resAAName << " is not a canonical amino acid type. Skipping contact" << endl;
            continue;
        }
        res_t aaType = SeqTools::aaToIdx(resAAName);
        Residue* targetRes = resPair.first;

        timer.start();
        residueFrame* query = defineQuery(resPair.first,resPair.second);
        frameTables[aaType]->setSearchParams(posCut,oriCut);
        allContactAADist.push_back(frameTables[aaType]->findBinderResAADist(query,targetRes->getResidueIndex()));
        timer.stop();

        mstreal duration = timer.getDuration(MstTimer::timeUnits::msec);
        cout << "Found distribution in " << duration << " ms" << endl;
        delete query;
    }
    return getJointProbability(allContactAADist);
}

vector<mstreal> binderScorer::getJointProbability(vector<vector<mstreal>> aaDists) {
    for (auto aaDist : aaDists) if (aaDist.size() != SeqToolsExtension::numAA()) MstUtils::error("amino acid distribution has wrong number of states","binderScorer::getJointProbability");
    vector<mstreal> jointDist;

    // find joint probability
    mstreal normConst;
    for (int i = 0; i < SeqToolsExtension::numAA(); i++) {
        mstreal prob = 1.0;
        for (int j = 0; j < aaDists.size(); j++) {
            prob *= aaDists[j][i];
        }
        jointDist.push_back(prob);
        normConst += prob;
    }

    // renormalize
    for (int i = 0; i < jointDist.size(); i++) jointDist[i] = jointDist[i]/normConst;
    return jointDist;
}