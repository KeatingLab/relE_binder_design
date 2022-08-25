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

binderScorer::binderScorer(const binderScorerParams& _params, Structure& _target) : params(_params), target(_target) {
    aaTypes = SeqToolsExtension::getAANames();
    targetName = MstSys::splitPath(target.getName(),1);
    pConts.load2DProbabilityDensities(params.potentialContactsJSONPath);
    setBackgroundSurfaceProbabilities();
}

binderScorer::binderScorer(const binderScorerParams& _params, Structure& complex, string binderChainIDsString, string targetChainIDsString) : params(_params) {
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
    if (targetResidues.empty()) MstUtils::error("No target residues found","binderScorer::binderScorer");
    target = Structure(targetResidues);
    targetName = MstSys::splitPath(complex.getName(),1);

    // extract binder chains from complex
    vector<Residue*> binderResidues;
    for (string binderChainID : binderChainIDs) {
        vector<Residue*> binderChainResidues = complex.getChainByID(binderChainID)->getResidues();
        binderResidues.insert(binderResidues.end(),binderChainResidues.begin(),binderChainResidues.end());
        // for (Residue* R : binderResidues) R->setName("UNK");
    }
    if (binderResidues.empty()) MstUtils::error("No binder residues found","binderScorer::binderScorer");
    binder = new Structure(binderResidues);
    binder->setName(targetName);
    complexMode = true; // will need to remember to delete this chain
    
    // define interface by looking at VDW contacts
    defineInterfaceUsingVDWContacts();

    // depending on what target residues are defined as the binding site, renormalize the frame probabilities
    // setFrameProbabilityTables();

    // in case the user will switch to using potential contacts later, load that data
    pConts.load2DProbabilityDensities(params.potentialContactsJSONPath);
    setBackgroundSurfaceProbabilities();
}

void binderScorer::setBinder(Structure* _binder) {
    binder = _binder;
}

void binderScorer::setTargetBindingSiteResidues(vector<Residue*> sel) {
    targetBindingResidues.clear();
    for (Residue* R : sel) targetBindingResidues.insert(R);
    cout << "Set " << targetBindingResidues.size() << " target residues" << endl;
    // setFrameProbabilityTables();
    pConts.setTargetResidues(sel);
}

void binderScorer::defineTargetBindingSiteResiduesByrSASA(mstreal relSASAthreshold) {
    targetBindingResidues.clear();
    sasaCalculator calc(target);
    bool relative = true;
    map<Residue*,mstreal> relSASA = calc.getResidueSASA(relative);
    for (auto it : relSASA) {
        if (it.second >= relSASAthreshold) targetBindingResidues.insert(it.first);
    }
    cout << "Set " << targetBindingResidues.size() << " target residues with relSASA threshold " << relSASAthreshold << endl;
    // setFrameProbabilityTables();
    vector<Residue*> targetResVec(targetBindingResidues.begin(),targetBindingResidues.end());
    pConts.setTargetResidues(targetResVec);
}

int binderScorer::defineInterfaceByPotentialContacts() {
    if (targetBindingResidues.empty()) MstUtils::error("Must define target residues before interface","binderScorer::defineInterfaceByPotentialContacts");
    if (binder->getResidues().empty()) MstUtils::error("Must define binder residues before interface","binderScorer::defineInterfaceByPotentialContacts");
    pConts.setBinderResidues(getBinderResidues());
    interfaceResidues = pConts.getContacts(false);
    cout << "Defined " << interfaceResidues.size() << " potential contacts" << endl;
    for (auto pair : interfaceResidues) {
        Residue* R = pair.first;
        // check that there are no target residues in the contacts that remain undefined
        if (targetBindingResidues.find(R) == targetBindingResidues.end()) MstUtils::error("Could not find target residue: "+R->getChainID()+MstUtils::toString(R->getNum()),"residueFrameBinderScorer::defineInterfaceByPotentialContacts");
        // check that there are no target residues in the contacts with an unrecognized amino acid type
        if (backgroundSurfaceProbabilities.find(SeqTools::aaToIdx(R->getName())) == backgroundSurfaceProbabilities.end()) MstUtils::error("Recognized residue types do not include target residue type: "+R->getName(),"binderScorer::defineInterfaceByPotentialContacts");
    }
    return interfaceResidues.size();
}

void binderScorer::defineInterfaceUsingVDWContacts() {
    if (!complexMode) MstUtils::error("Can only define interface using van der Waals contacts when scoring a full-atom complex","binderScorer::defineInterfaceUsingVDWContacts");
    vdwContacts C(target.getResidues(),binder->getResidues());
    interfaceResidues = C.getInteractingResPairs();
    for (auto pair : interfaceResidues) targetBindingResidues.insert(pair.first);

    cout << interfaceResidues.size() << " interface contacts when defined by van der Waals interactions" << endl;
}

void binderScorer::setBackgroundSurfaceProbabilities() {
    backgroundSurfaceProbabilities[SeqTools::aaToIdx("ALA")] = 0.048049;
    backgroundSurfaceProbabilities[SeqTools::aaToIdx("ARG")] = 0.090171;
    backgroundSurfaceProbabilities[SeqTools::aaToIdx("ASN")] = 0.062437;
    backgroundSurfaceProbabilities[SeqTools::aaToIdx("ASP")] = 0.095417;
    backgroundSurfaceProbabilities[SeqTools::aaToIdx("CYS")] = 0.002533;
    backgroundSurfaceProbabilities[SeqTools::aaToIdx("GLN")] = 0.070898;
    backgroundSurfaceProbabilities[SeqTools::aaToIdx("GLU")] = 0.125364;
    backgroundSurfaceProbabilities[SeqTools::aaToIdx("GLY")] = 0.029287;
    backgroundSurfaceProbabilities[SeqTools::aaToIdx("HIS")] = 0.023232;
    backgroundSurfaceProbabilities[SeqTools::aaToIdx("ILE")] = 0.016676;
    backgroundSurfaceProbabilities[SeqTools::aaToIdx("LEU")] = 0.033629;
    backgroundSurfaceProbabilities[SeqTools::aaToIdx("LYS")] = 0.132090;
    backgroundSurfaceProbabilities[SeqTools::aaToIdx("MET")] = 0.009195;
    backgroundSurfaceProbabilities[SeqTools::aaToIdx("PHE")] = 0.048049;
    backgroundSurfaceProbabilities[SeqTools::aaToIdx("PRO")] = 0.012781;
    backgroundSurfaceProbabilities[SeqTools::aaToIdx("SER")] = 0.068929;
    backgroundSurfaceProbabilities[SeqTools::aaToIdx("THR")] = 0.058319;
    backgroundSurfaceProbabilities[SeqTools::aaToIdx("TRP")] = 0.006683;
    backgroundSurfaceProbabilities[SeqTools::aaToIdx("TYR")] = 0.016602;
    backgroundSurfaceProbabilities[SeqTools::aaToIdx("VAL")] = 0.027414;
}

/* --- --- --- --- --- residueBackboneBinderScorer --- --- --- --- --- */

mstreal residueBackboneBinderScorer::scoreBinder() {
    binderScore = 0.0;
    binderScorePerContact.clear();
    numMatchesPerContact.clear();
    numNativeMatchesPerContact.clear();
    nMatchesAAPerContact.clear();
    searchTimePerContact.clear();
    backboneAtomDistancesPerContact.clear();

    if (interfaceResidues.empty()) {
        cout << "No contacts to score." << endl;
        binderScore = -1000.0;
        return binderScore;
    }

    // score each pair of residues with potential to contact
    for (auto contactingResidues : interfaceResidues) {
        int totalMatches = 0;
        int nativeAAMatches = 0;
        mstreal contactScore = 0.0;
        scoreContact(contactingResidues,totalMatches,nativeAAMatches,contactScore);
        binderScore+=contactScore;
    }
    return binderScore;
}

void residueBackboneBinderScorer::scoreContact(pair<Residue*,Residue*> contactingRes, int& totalMatches, int& nativeMatches, mstreal& contactScore) {
    resPairSearcher.setQuery(contactingRes.first,contactingRes.second);
    backboneAtomDistancesPerContact.push_back(resPairSearcher.getQuery().getAllbbAtomDistances());
    timer.start();
    totalMatches = resPairSearcher.searchForMatches();
    timer.stop();
    mstreal searchTime = timer.getDuration(MstTimer::timeUnits::msec);
    cout << "Took " << searchTime << " (ms) to find " << totalMatches << " matches." << endl;
    if (totalMatches > 0) { // the exact threshold should be determined more rigorously
        res_t targetResType = SeqTools::aaToIdx(contactingRes.first->getName());
        nativeMatches = resPairSearcher.getNumMatchesWithResidueType(true);
        mstreal nativeResProbFromMatches = mstreal(nativeMatches+1)/mstreal(totalMatches+aaTypes.size());
        contactScore = -log2(nativeResProbFromMatches/backgroundSurfaceProbabilities[targetResType]);
        nMatchesAAPerContact.push_back(resPairSearcher.getNumMatchesByAAType(true));
    } else {
        nativeMatches = 0;
        contactScore = 10.0;
        nMatchesAAPerContact.push_back(vector<int>(20,0));
    }
    binderScorePerContact.push_back(contactScore);
    numMatchesPerContact.push_back(totalMatches);
    numNativeMatchesPerContact.push_back(nativeMatches);
    searchTimePerContact.push_back(searchTime);
}

void residueBackboneBinderScorer::writeBinderScoresToFile(bool append) {
    if (binder_info_out == nullptr) {
        binder_info_out = new fstream;
        string filename = targetName + "_binderScores.tsv";
        MstUtils::openFile(*binder_info_out,filename,fstream::out);
        *binder_info_out << "targetName\tbinderName\tbinderScore\tnumPotentialContacts\tnumBinderRes" << endl;
    }
    *binder_info_out << targetName << "\t";
    *binder_info_out << binder->getName() << "\t";
    *binder_info_out << binderScore << "\t";
    *binder_info_out << countDesignableContacts() << "\t";
    *binder_info_out << binder->residueSize();
    *binder_info_out << endl;
}

void residueBackboneBinderScorer::writeContactScoresToFile(bool append) {
    if (contact_info_out == nullptr) {
        contact_info_out = new fstream;
        string filename = targetName + "_contactScores.tsv";
        MstUtils::openFile(*contact_info_out,filename,fstream::out);
        *contact_info_out << "targetName\tbinderName\ttargetRes\ttargetResName\tbinderRes\tbinderResName\tbinderScore\tnTotalMatches\tnNativeMatches\tdistanceCa\tsearchTime_ms" << endl;
    }
    for (int i = 0; i < interfaceResidues.size(); i++) {
        Residue* targetR = interfaceResidues[i].first;
        Residue* binderR = interfaceResidues[i].second;
        mstreal distanceCA = (targetR->findAtom("CA")->getCoor() - binderR->findAtom("CA")->getCoor()).norm();
        *contact_info_out << targetName << "\t";
        *contact_info_out << binderR->getStructure()->getName() << "\t";
        *contact_info_out << targetR->getChainID() << targetR->getNum() << "\t" << targetR->getName() << "\t";
        *contact_info_out << binderR->getChainID() << binderR->getNum() << "\t" << binderR->getName()<< "\t";
        *contact_info_out << binderScorePerContact[i] << "\t";
        *contact_info_out << numMatchesPerContact[i] << "\t";
        *contact_info_out << numNativeMatchesPerContact[i] << "\t";
        *contact_info_out << distanceCA << "\t";
        *contact_info_out << searchTimePerContact[i];
        *contact_info_out << endl;
    }
}

void residueBackboneBinderScorer::writeTrainingDataToFile(bool append) {
    if (training_info_out == nullptr) {
        training_info_out = new fstream;
        string filename = targetName + "_trainingData.csv";
        MstUtils::openFile(*training_info_out,filename,fstream::out);
        // complex and contacting residues
        *training_info_out << "targetName,binderName,targetRes,targetResName,binderRes,binderResName,";

        // geometric features describing contacting residues
        *training_info_out << "N-N,N-Ca,N-C,N-O,Ca-N,Ca-Ca,Ca-C,Ca-O,C-N,C-Ca,C-C,C-O,O-N,O-Ca,O-C,O-O,";

        // statistics from the database
        *training_info_out << "binderScore,nTotalMatches,nNativeMatches,";

        // counts of each aa type
        *training_info_out << "A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y" << endl;
    }
    for (int i = 0; i < interfaceResidues.size(); i++) {
        // complex and contacting residues
        Residue* targetR = interfaceResidues[i].first;
        Residue* binderR = interfaceResidues[i].second;
        mstreal distanceCA = (targetR->findAtom("CA")->getCoor() - binderR->findAtom("CA")->getCoor()).norm();
        *training_info_out << targetName << ",";
        *training_info_out << binderR->getStructure()->getName() << ",";
        *training_info_out << targetR->getChainID() << targetR->getNum() << "," << targetR->getName() << ",";
        *training_info_out << binderR->getChainID() << binderR->getNum() << "," << binderR->getName()<< ",";

        // geometric features describing contacting residues
        CartesianPoint backboneAtomDistances = backboneAtomDistancesPerContact[i];
        for (mstreal distance : backboneAtomDistances) *training_info_out << distance << ",";

        // statistics from the database
        *training_info_out << binderScorePerContact[i] << ",";
        *training_info_out << numMatchesPerContact[i] << ",";
        *training_info_out << numNativeMatchesPerContact[i];

        // counts of each aa type
        for (int count : nMatchesAAPerContact[i]) *training_info_out << "," << count;
        *training_info_out << endl;
    }
}


/* --- --- --- --- --- residueFrameBinderScorer --- --- --- --- --- */
// residueFrameBinderScorer::residueFrameBinderScorer(const binderScorerParams& params, augmentedStructure& _target) : binderScorer(params,_target), augmentedTarget(_target) {
//     prepareVoxelGrids(params.frameDBPath);
// }

// residueFrameBinderScorer::residueFrameBinderScorer(const binderScorerParams& params, augmentedStructure& complex, string binderChainIDsString, string targetChainIDsString) : binderScorer(params,complex,binderChainIDsString,targetChainIDsString) {
//     augmentedBinder = new augmentedStructure(binder->getResidues());
//     augmentedBinder->setName(targetName);
    
//     prepareVoxelGrids(params.frameDBPath);

//     // define interface by looking at VDW contacts
//     defineInterfaceUsingVDWContacts();

//     // in case the user will switch to using potential contacts later, load that data
//     pConts.load2DProbabilityDensities(params.potentialContactsJSONPath);
//     setBackgroundSurfaceProbabilities();
// }

// void residueFrameBinderScorer::setBinder(augmentedStructure* _binder) {
//     if (complexMode) MstUtils::error("Cannot set another binder in complex mode","residueFrameBinderScorer::setBinder");
//     binder = _binder;
//     // for (Residue* R : binder->getResidues()) R->setName("UNK");
//     pConts.setBinderResidues(binder->getResidues());

//     // reset all variables
//     designabilityScore = 0.0;
//     numNonDesignable = 0;
// }

// void residueFrameBinderScorer::setTargetBindingSiteResidues(vector<Residue*> sel) {
//     targetBindingResidues.clear();
//     for (Residue* R : sel) targetBindingResidues.insert(R);
//     cout << "Set " << targetBindingResidues.size() << " target residues" << endl;
//     // setFrameProbabilityTables();
//     pConts.setTargetResidues(sel);
// }

// void residueFrameBinderScorer::defineTargetBindingSiteResiduesByrSASA(mstreal relSASAthreshold) {
//     targetBindingResidues.clear();
//     sasaCalculator calc(target);
//     bool relative = true;
//     map<Residue*,mstreal> relSASA = calc.getResidueSASA(relative);
//     for (auto it : relSASA) {
//         if (it.second >= relSASAthreshold) targetBindingResidues.insert(it.first);
//     }
//     cout << "Set " << targetBindingResidues.size() << " target residues with relSASA threshold " << relSASAthreshold << endl;
//     // setFrameProbabilityTables();
//     vector<Residue*> targetResVec(targetBindingResidues.begin(),targetBindingResidues.end());
//     pConts.setTargetResidues(targetResVec);
// }

// int residueFrameBinderScorer::defineInterfaceByPotentialContacts() {
//     if (targetBindingResidues.empty()) MstUtils::error("Must define target residues before interface","residueFrameBinderScorer::defineInterfaceByPotentialContacts");
//     pConts.setBinderResidues(getBinderResidues());
//     interfaceResidues = pConts.getContacts(false);
//     cout << "Defined " << interfaceResidues.size() << " potential contacts" << endl;
//     // check that there are no target residues in the contacts that remain undefined
//     for (auto pair : interfaceResidues) {
//         Residue* R = pair.first;
//         if (targetBindingResidues.find(R) == targetBindingResidues.end()) MstUtils::error("Could not find target residue: "+R->getChainID()+MstUtils::toString(R->getNum()),"residueFrameBinderScorer::defineInterfaceByPotentialContacts");
//     }
//     interfaceCounts.clear();
//     interfaceDesignabilityScores.clear();
//     interfaceSequenceCompatibilityScores.clear();
//     interfaceTotalNumberOfMatches.clear();
//     searchTimes.clear();
//     return interfaceResidues.size();
// }

// mstreal residueFrameBinderScorer::scoreInterfaceDesignability() {
//     if (!interfaceDesignabilityScores.empty()) return designabilityScore;
//     else if (interfaceResidues.empty()) {
//         cout << "No contacts to score." << endl;
//         return 0.0;
//     }

//     // score each pair of residues with potential to contact
//     for (auto resPair : interfaceResidues) {
//         int residueFramePairCounts = residueFrameCounts(resPair);
//         mstreal residuePairScore = logCounts(residueFramePairCounts);
//         designabilityScore += residuePairScore;
//         interfaceCounts.push_back(residueFramePairCounts);
//         interfaceDesignabilityScores.push_back(residuePairScore);
//     }
//     nonDesignableInterfaceResidues = pConts.getNonDesignableContacts();
//     numNonDesignable = nonDesignableInterfaceResidues.size();
//     return designabilityScore + nonDesignableContactPenalty*numNonDesignable;
// }

// mstreal residueFrameBinderScorer::scoreInterfaceSequenceCompatibility() {
//     if (!interfaceSequenceCompatibilityScores.empty()) return sequenceCompatibilityScore;
//     else if (interfaceResidues.empty()) {
//         cout << "No contacts to score." << endl;
//         return 0.0;
//     }

//     // score each pair of residues with potential to contact
//     for (auto resPair : interfaceResidues) {
//         int totalMatches = 0;
//         int nativeAAMatches = 0;
//         mstreal sequenceCompatibilityContactScore = 0.0;
//         residueFrameSequenceCompatibility(resPair,totalMatches,nativeAAMatches,sequenceCompatibilityContactScore);
//         sequenceCompatibilityScore+=sequenceCompatibilityContactScore;
//         interfaceSequenceCompatibilityScores.push_back(sequenceCompatibilityContactScore);
//         interfaceTotalNumberOfMatches.push_back(totalMatches - 20);
//     }
//     nonDesignableInterfaceResidues = pConts.getNonDesignableContacts();
//     numNonDesignable = nonDesignableInterfaceResidues.size();
//     return sequenceCompatibilityScore; // + nonDesignableContactPenalty*numNonDesignable;
// }

// void residueFrameBinderScorer::writeBinderDesignabilityScoresToFile(bool append) {
//     if (binder_info_out == nullptr) {
//         binder_info_out = new fstream;
//         string filename = targetName + "_binderScores.tsv";
//         MstUtils::openFile(*binder_info_out,filename,fstream::out);
//         *binder_info_out << "targetName\tbinderName\tbinderScore\tresidueFrameScore\tnumPotentialContacts\tnumNonDesignableContacts\tnumBinderRes" << endl;
//     }
//     *binder_info_out << targetName << "\t";
//     *binder_info_out << binder->getName() << "\t";
//     *binder_info_out << designabilityScore + nonDesignableContactPenalty*numNonDesignable << "\t";
//     *binder_info_out << designabilityScore << "\t";
//     *binder_info_out << countDesignableContacts() << "\t";
//     *binder_info_out << numNonDesignable << "\t";
//     *binder_info_out << binder->residueSize();
//     *binder_info_out << endl;
// }

// void residueFrameBinderScorer::writeContactDesignabilityScoresToFile(bool append) {
//     if (contact_info_out == nullptr) {
//         contact_info_out = new fstream;
//         string filename = targetName + "_contactScores.tsv";
//         MstUtils::openFile(*contact_info_out,filename,fstream::out);
//         *contact_info_out << "targetName\tbinderName\ttargetRes\ttargetResName\tbinderRes\tbinderResName\tresidueFrameCounts\tresidueFrameScore\tdistanceCa\tsearchTime_ms" << endl;
//     }
//     for (int i = 0; i < interfaceResidues.size(); i++) {
//         Residue* targetR = interfaceResidues[i].first;
//         Residue* binderR = interfaceResidues[i].second;
//         mstreal distanceCA = (targetR->findAtom("CA")->getCoor() - binderR->findAtom("CA")->getCoor()).norm();
//         *contact_info_out << targetName << "\t";
//         *contact_info_out << binderR->getStructure()->getName() << "\t";
//         *contact_info_out << targetR->getChainID() << targetR->getNum() << "\t" << targetR->getName() << "\t";
//         *contact_info_out << binderR->getChainID() << binderR->getNum() << "\t" << binderR->getName()<< "\t";
//         *contact_info_out << interfaceCounts[i] << "\t";
//         *contact_info_out << interfaceDesignabilityScores[i] << "\t";
//         *contact_info_out << distanceCA << "\t";
//         *contact_info_out << searchTimes[i];
//         *contact_info_out << endl;
//     }
// }

// void residueFrameBinderScorer::writeBinderSequenceCompatibilityScoresToFile(bool append) {
//     if (binder_info_out == nullptr) {
//         binder_info_out = new fstream;
//         string filename = targetName + "_binderScores.tsv";
//         MstUtils::openFile(*binder_info_out,filename,fstream::out);
//         *binder_info_out << "targetName\tbinderName\tsequenceCompatiblityScore\tnumPotentialContacts\tnumBinderRes" << endl;
//     }
//     *binder_info_out << targetName << "\t";
//     *binder_info_out << binder->getName() << "\t";
//     *binder_info_out << sequenceCompatibilityScore << "\t";
//     *binder_info_out << countDesignableContacts() << "\t";
//     *binder_info_out << binder->residueSize();
//     *binder_info_out << endl;
// }

// void residueFrameBinderScorer::writeContactSequenceCompatibilityScoresToFile(bool append) {
//     if (contact_info_out == nullptr) {
//         contact_info_out = new fstream;
//         string filename = targetName + "_contactScores.tsv";
//         MstUtils::openFile(*contact_info_out,filename,fstream::out);
//         *contact_info_out << "targetName\tbinderName\ttargetRes\ttargetResName\tbinderRes\tbinderResName\tsequenceCompatibilityScore\tnTotalMatches\tdistanceCa\tsearchTime_ms" << endl;
//     }
//     for (int i = 0; i < interfaceResidues.size(); i++) {
//         Residue* targetR = interfaceResidues[i].first;
//         Residue* binderR = interfaceResidues[i].second;
//         mstreal distanceCA = (targetR->findAtom("CA")->getCoor() - binderR->findAtom("CA")->getCoor()).norm();
//         *contact_info_out << targetName << "\t";
//         *contact_info_out << binderR->getStructure()->getName() << "\t";
//         *contact_info_out << targetR->getChainID() << targetR->getNum() << "\t" << targetR->getName() << "\t";
//         *contact_info_out << binderR->getChainID() << binderR->getNum() << "\t" << binderR->getName()<< "\t";
//         *contact_info_out << interfaceSequenceCompatibilityScores[i] << "\t";
//         *contact_info_out << interfaceTotalNumberOfMatches[i] << "\t";
//         *contact_info_out << distanceCA << "\t";
//         *contact_info_out << searchTimes[i];
//         *contact_info_out << endl;
//     }
// }

// void residueFrameBinderScorer::writeContactPropertyToFile(string dirPath, bool designabilityScore) {
//     if (designabilityScore && interfaceDesignabilityScores.empty()) MstUtils::error("Cannot write designability scores before computing","residueFrameBinderScorer::writeContactPropertyToFile");
//     if (!designabilityScore && interfaceSequenceCompatibilityScores.empty()) MstUtils::error("Cannot write sequence compatibility scores before computing","residueFrameBinderScorer::writeContactPropertyToFile");
//     fstream conts_out;
//     if (dirPath != "") dirPath = dirPath + "/";
//     string filename = dirPath + binder->getName() + "_contactScoresSimple.tsv";
//     cout << "filename: " << filename << endl;
//     MstUtils::openFile(conts_out, filename, fstream::out, "residueFrameBinderScorer::writeContactPropertyToFile");

//     for (int i = 0; i < interfaceResidues.size(); i++) {
//         Residue* Ri = interfaceResidues[i].first;
//         Residue* Rj = interfaceResidues[i].second;
//         mstreal val;
//         if (designabilityScore) val = interfaceDesignabilityScores[i];
//         else val = interfaceSequenceCompatibilityScores[i];
//         conts_out << Ri->getChainID() << Ri->getNum() << "\t";
//         conts_out << Rj->getChainID() << Rj->getNum() << "\t";
//         conts_out << val;
//         conts_out << endl;
//     }
// }

// void residueFrameBinderScorer::writeResidueClashesToFile(bool append) {
//     if (clash_info_out == nullptr) {
//         clash_info_out = new fstream;
//         string filename = targetName + "_residueClashes.tsv";
//         MstUtils::openFile(*clash_info_out,filename,fstream::out);
//         *clash_info_out << "targetName\tbinderName\ttargetRes\ttargetResName\tbinderRes\tbinderResName\tdistance\tnormCbDistance\tsearchTime_ms" << endl;
//     }
//     vector<pair<Residue*,Residue*>> residueClashes = pConts.getNonDesignableContacts();
//     for (int i = 0; i < residueClashes.size(); i++) {
//         Residue* targetR = residueClashes[i].first;
//         Residue* binderR = residueClashes[i].second;
//         mstreal distanceCA = (targetR->findAtom("CA")->getCoor() - binderR->findAtom("CA")->getCoor()).norm();
//         mstreal normCbDist = pConts.getNormalizedCbDistance(targetR,binderR);
//         *clash_info_out << targetName << "\t";
//         *clash_info_out << binderR->getStructure()->getName() << "\t";
//         *clash_info_out << targetR->getChainID() << targetR->getNum() << "\t" << targetR->getName() << "\t";
//         *clash_info_out << binderR->getChainID() << binderR->getNum() << "\t" << binderR->getName()<< "\t";
//         *clash_info_out << distanceCA << "\t" << normCbDist;
//         *clash_info_out << endl;
//     }
// }

// void residueFrameBinderScorer::writeResiduesClashesToSimpleFile(string dirPath) {
//     fstream conts_out;
//     if (dirPath != "") dirPath = dirPath + "/";
//     string filename = dirPath + MstSys::splitPath(binder->getName(),1) + "_residueClashesSimple.tsv";
//     cout << "filename: " << filename << endl;
//     MstUtils::openFile(conts_out, filename, fstream::out, "residueFrameBinderScorer::writeContactPropertyToFile");

//     for (int i = 0; i < interfaceResidues.size(); i++) {
//         Residue* Ri = interfaceResidues[i].first;
//         Residue* Rj = interfaceResidues[i].second;
//         conts_out << Ri->getChainID() << Ri->getNum() << "\t";
//         conts_out << Rj->getChainID() << Rj->getNum() << "\t";
//         conts_out << 1.0;
//         conts_out << endl;
//     }
// }

// void residueFrameBinderScorer::writePSSM(string dirPath) {
//     fstream pssm_out;
//     if (dirPath != "") dirPath = dirPath + "/";
//     string filename = dirPath + binder->getName() + "_binderPSSM.csv";
//     cout << "filename: " << filename << endl;
//     MstUtils::openFile(pssm_out, filename, fstream::out, "interfaceScorer::writePSSM");
    
//     // write header
//     set<string> aaNames = SeqToolsExtension::getAANames();
//     vector<string> aaNamesVec(aaNames.begin(),aaNames.end());
//     pssm_out << "position," << MstUtils::join(",",aaNamesVec) << endl;    

//     for (Residue* binderR : binder->getResidues()) {
//         // if residue has potential contacts, we can estimate the amino acid probability distribution
//         vector<mstreal> aaDist = getAADistAtPos(binderR);
//         if (aaDist.empty()) continue;
//         pssm_out << binderR->getResidueIndexInChain()+1; // seeds are 1-indexed
//         for (mstreal val : aaDist) {
//             pssm_out << "," << val;
//         }
//         pssm_out << endl;
//     }
//     pssm_out.close();
// }

// void residueFrameBinderScorer::prepareVoxelGrids(string frameDBPath) {
//     timer.start();

//     frameDB DB(frameDBPath);
//     vector<mobileFrame*> loaded = DB.loadAllFrames();
//     cout << "Loaded " << loaded.size() << " frames in total." << endl;

//     for (string aa : aaTypes) {
//         cout << "Loading frames with amino acid type: " << aa << endl;
//         vector<mobileFrame*> currentFrames;
//         res_t current_aa = SeqTools::aaToIdx(aa);
//         for (mobileFrame* frame : loaded) {
//             if (frame->getResIIndex() == current_aa) {
//                 currentFrames.push_back(frame);
//             }
//         }
//         boundingBox bbox(1.0);
//         bbox.update(currentFrames);
//         bbox.printBounds();
//         frameTable* mobileResFrames = new frameTable(bbox, posCut/mstreal(2), 360/int(oriCut/mstreal(2)),verbose);
//         for (mobileFrame* frame : currentFrames) mobileResFrames->insertFrame(frame);
//         frameTables[current_aa] = mobileResFrames;
//     }

//     timer.stop();
//     cout << "It took " << timer.getDuration() << " seconds to load all frames into hash tables" << endl;
// }

// // void binderScorer::setFrameProbabilityTables() {
// //     // for each target residue with the same amino acid identity, compute the occupied volume
// //     for (string aa : aaTypes) {
// //         cout << "Set target residues with identity: " << aa << endl;
// //         res_t current_aa = SeqTools::aaToIdx(aa);
// //         frameProbability* fTable = frameTables.at(current_aa);
// //         for (Residue* targetRes : targetResidues) {
// //             if (SeqTools::aaToIdx(targetRes->getName()) == current_aa) {
// //                 if (fTable->isTargetResidueDefined(targetRes->getResidueIndex())) continue;
// //                 fTable->setTargetResidue(targetRes->getResidueIndex(),targetBackbone,renormalizeProbabilities);
// //             }
// //         }
// //     }
// // }

// int residueFrameBinderScorer::residueFrameCounts(pair<Residue*,Residue*> resPair) {
//     string resAAName = resPair.first->getName();
//     if (aaTypes.find(resAAName) == aaTypes.end()) {
//         cout << "Warning: amino acid type: " << resAAName << " is not a canonical amino acid type. Skipping contact" << endl;
//         return 0.0;
//     }
//     res_t aaType = SeqTools::aaToIdx(resAAName);
//     Residue* targetRes = resPair.first;
//     Residue* binderRes = resPair.second;
//     if (verbose) cout << "Searching contact between " << targetRes->getChainID() << targetRes->getNum() << " and " << binderRes->getChainID() << binderRes->getNum() << endl;
    
//     timer.start();
//     residueFrame* query = defineQuery(targetRes,binderRes);
//     // frameTables[aaType]->setSearchParams(posCut,oriCut);
//     int counts = frameTables[aaType]->countSimilarFrames(query,posCut,oriCut);
//     timer.stop();

//     mstreal duration = timer.getDuration(MstTimer::timeUnits::msec);
//     if (verbose) cout << "Found " << counts << " matches." << endl;
//     searchTimes.push_back(duration);

//     delete query;
//     return counts;
// }

// void residueFrameBinderScorer::residueFrameSequenceCompatibility(pair<Residue*,Residue*> resPair, int &totalMatches, int &nativeAAMatches, mstreal &sequenceCompatibilityScore) {
//     string resAAName = resPair.first->getName();
//     if (aaTypes.find(resAAName) == aaTypes.end()) {
//         cout << "Warning: amino acid type: " << resAAName << " is not a canonical amino acid type. Skipping contact" << endl;
//         totalMatches = 0;
//         nativeAAMatches = 0;
//         sequenceCompatibilityScore = 0.0;
//     }
//     res_t aaType = SeqTools::aaToIdx(resAAName);
//     Residue* targetRes = resPair.first;
//     Residue* binderRes = resPair.second;
//     if (verbose) cout << "Searching contact between " << targetRes->getChainID() << targetRes->getNum() << " and " << binderRes->getChainID() << binderRes->getNum() << endl;
    
//     timer.start();
//     residueFrame* query = defineQuery(targetRes,binderRes);
//     for (auto const& table : frameTables) {
//         int counts = frameTables[table.first]->countSimilarFrames(query,posCut,oriCut);
//         totalMatches += counts;
//         totalMatches+=1; //psuedocount
//         if (table.first == aaType) {
//             nativeAAMatches += counts;
//             nativeAAMatches+=1; //pseudocount
//         }
//     }
//     timer.stop();

//     // Calculate score
//     cout << "totalMatches " << totalMatches << endl;
//     cout << "nativeAAMatches " << nativeAAMatches << endl;
//     cout << "bgProb " << backgroundSurfaceProbabilities[aaType] << endl;
//     sequenceCompatibilityScore = -log2((mstreal(nativeAAMatches)/mstreal(totalMatches))/backgroundSurfaceProbabilities[aaType]);

//     mstreal duration = timer.getDuration(MstTimer::timeUnits::msec);
//     if (verbose) cout << "Found " << totalMatches << " total matches and " << nativeAAMatches << " with the native amino acid" << endl;
//     searchTimes.push_back(duration);

//     delete query;
// }

// mstreal residueFrameBinderScorer::logProb(mstreal numerator, mstreal denominator, mstreal pseudocount) {
//     return log(numerator+pseudocount) - log(denominator+pseudocount);
// }

// vector<mstreal> residueFrameBinderScorer::getAADistAtPos(Residue* binderRes) {
//     vector<mstreal> aaDist;

//     // collect the interacting pairs
//     vector<pair<Residue*,Residue*>> contactsWithBinderR;
//     for (auto pair : interfaceResidues) {
//         if (pair.second == binderRes) {
//             contactsWithBinderR.push_back(pair);
//         }
//     }
//     if (contactsWithBinderR.empty()) return aaDist;

//     // estimate the 'diagonal' of the joint distribution 
//     // e.g. what is the probability of some amino acid at this binder residue, given
//     // each of the target residues (assuming independence)
//     vector<vector<mstreal>> allContactAADist;
//     for (pair<Residue*,Residue*> resPair : contactsWithBinderR) {
//         string resAAName = resPair.first->getName();
//         if (aaTypes.find(resAAName) == aaTypes.end()) {
//             cout << "Warning: amino acid type: " << resAAName << " is not a canonical amino acid type. Skipping contact" << endl;
//             continue;
//         }
//         res_t aaType = SeqTools::aaToIdx(resAAName);
//         Residue* targetRes = resPair.first;

//         timer.start();
//         residueFrame* query = defineQuery(resPair.first,resPair.second);
//         allContactAADist.push_back(frameTables[aaType]->findBinderResAADist(query,posCut,oriCut));
//         timer.stop();

//         mstreal duration = timer.getDuration(MstTimer::timeUnits::msec);
//         cout << "Found distribution in " << duration << " ms" << endl;
//         delete query;
//     }
//     return getJointProbability(allContactAADist);
// }

// vector<mstreal> residueFrameBinderScorer::getJointProbability(vector<vector<mstreal>> aaDists) {
//     for (auto aaDist : aaDists) if (aaDist.size() != SeqToolsExtension::numAA()) MstUtils::error("amino acid distribution has wrong number of states","binderScorer::getJointProbability");
//     vector<mstreal> jointDist;

//     // find joint probability
//     mstreal normConst;
//     for (int i = 0; i < SeqToolsExtension::numAA(); i++) {
//         mstreal prob = 1.0;
//         for (int j = 0; j < aaDists.size(); j++) {
//             prob *= aaDists[j][i];
//         }
//         jointDist.push_back(prob);
//         normConst += prob;
//     }

//     // renormalize
//     for (int i = 0; i < jointDist.size(); i++) jointDist[i] = jointDist[i]/normConst;
//     return jointDist;
// }

// /* --- --- --- --- --- binderBackboneScorer --- --- --- --- --- */

// mstreal binderBackboneScorer::scoreBackbone(Structure* S) {
//     if (!info_out.is_open()) openFile();
//     vector<Residue*> Sresidues = S->getResidues();
//     int segmentLength = segmentSearcher.getSegLen();
//     mstreal total_score = 0;
//     if (S->residueSize() < segmentLength) MstUtils::error("The length of segments in the library is longer than the structure that is to be scored","binderBackboneScorer::scoreBackbone");
//     for (int i = 0; i < S->residueSize() - segmentLength + 1; i++) {
//         Structure segment(vector<Residue*>(Sresidues.begin()+i,Sresidues.begin()+i+segmentLength));
//         mstreal RMSD = segmentSearcher.findLowestRMSDSegment(&segment);
//         int nMatches = segmentSearcher.getNumMatchesInDB();
//         mstreal score = -log(nMatches+1);
//         cout << "Found segment matching position " << i << " with RMSD: " << RMSD;
//         cout << " and nMatches: " << nMatches << endl;
//         info_out << name << "," << i << ",";
//         info_out << segmentLength << "," << nMatches << ",";
//         info_out << score << endl;
//     }
//     return total_score;
// };

// void binderBackboneScorer::openFile() {
//     string fileName = name + "_backboneSegmentScore.csv";
//     MstUtils::openFile(info_out,fileName,fstream::out);
//     info_out << "name,nTermResIdx,segLen,nMatches,score" << endl;
// };

// /* --- --- --- --- --- binderBackboneScorer --- --- --- --- --- */