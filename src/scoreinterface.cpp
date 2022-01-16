#include "scoreinterface.h"

/* --- --- --- --- --- interfaceScorer --- --- --- --- --- */

interfaceScorer::interfaceScorer(string frameDBPath, mstreal _posCut, mstreal _oriCut) {
    if (_posCut > 0) posCut = _posCut;
    if (_oriCut > 0) oriCut = _oriCut;
    setParams();
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
        vector<mstreal> bbox = {0,0,0,0,0,0};
        bool verbose = true;
        positionHasher::updateBoundingBox(currentFrames,bbox,verbose);
        frameTable* fTable = new frameTable(bbox, 0.25, 36);
        for (mobileFrame* frame : currentFrames) fTable->insertFrame(frame);
        frameTables[current_aa] = fTable;
    }

    timer.stop();
    cout << "It took " << timer.getDuration() << " seconds to load all frames into hash tables" << endl;
}

void interfaceScorer::scoreInterface(vector<pair<Residue*,Residue*>> _interfaceResidues, alignInteractingFrames* alignF) {
    // reset all variables
    matchingFrames.clear();
    searchTimes.clear();

    // define interface, if necessary
    if (!_interfaceResidues.empty()) {
        interfaceResidues = _interfaceResidues;
    } else if (interfaceResidues.empty()) {
        cout << "No contacts provided, defining interface using simple distance check" << endl;
        defineInterface();
    }
    cout << interfaceResidues.size() << " interface contacts total" << endl;

    // score each interface residue pair
    for (auto resPair : interfaceResidues) {
        scoreResiduePair(resPair,alignF);
    }
}

void interfaceScorer::scoreResiduePair(pair<Residue*,Residue*> resPair, alignInteractingFrames* alignF) {
    string resAAName = resPair.first->getName();
    if (aaTypes.find(resAAName) == aaTypes.end()) {
        // MstUtils::error("Amino acid type:"+resAAName+" not a canonical amino acid");
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

void interfaceScorer::writeContactScoresToFile(string name, bool append) {
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
        *contact_info_out << complex->getName() << "\t";
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

void interfaceScorer::writeContactMatchesToFile(string name, bool append) {
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
            *match_info_out << complex->getName() << "\t";
            *match_info_out << targetR->getChainID() << targetR->getNum() << "\t" << targetR->getName() << "\t";
            *match_info_out << binderR->getChainID() << binderR->getNum() << "\t" << binderR->getName()<< "\t";
            *match_info_out << frame->getTarget() << "\t" << frame->getResI() << "\t" << frame->getResJ() << "\t";
            *match_info_out << SeqTools::idxToTriple(frame->getResJIndex());
            *match_info_out << endl;
        }
    }
}

void interfaceScorer::writeMatchStructures(string name, alignInteractingFrames& alignF) {
    fstream pdb_out;
    MstUtils::openFile(pdb_out, name, fstream::out, "interfaceScorer::writeMatchStructures");

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

void interfaceScorer::writeContactPropertyToFile(string name, interfaceScorer::property toWrite) {
    fstream conts_out;
    string filename = name + "_contactProperty.tsv";
    MstUtils::openFile(conts_out, filename, fstream::out, "interfaceScorer::writeContactPropertyToFile");

    conts_out << "RiChain\tRiNum\tRjChain\tRjNum\tVal" << endl;
    for (int i = 0; i < interfaceResidues.size(); i++) {
        Residue* Ri = interfaceResidues[i].first;
        Residue* Rj = interfaceResidues[i].second;
        conts_out << Ri->getChainID() << "\t" << Ri->getNum() << "\t";
        conts_out << Rj->getChainID() << "\t" << Rj->getNum() << "\t";
        mstreal val = 0.0;
        if (property::COVERAGE == toWrite) {
            val = (matchingFrames[i].empty()) ? 1 : -1; //negative score is considered more favorable, so -1 means covered
        } else {
            MstUtils::error("property type to write to file is not supported","interfaceScorer::writeContactPropertyToFile");
        }
        conts_out << val;
        conts_out << endl;
    }
}

void interfaceScorer::setParams() {
    aaTypes = {"ALA","ARG","ASN","ASP",
                "CYS","GLN","GLU","GLY",
                "HIS","ILE","LEU","LYS",
                "MET","PHE","PRO","SER",
                "THR","TRP","TYR","VAL"};
}

void interfaceScorer::defineInterface() {
    interfaceResidues.clear();
    mstreal cutoff = 7.5;

    // use a simple distance-based definition for now
    vector<Atom*> target, binder;
    for (Chain* C : proteinChains) {
        C->getResidues();
        for (Residue* R : C->getResidues()) {
            for (Atom* A : R->getAtoms()) {
                target.push_back(A);
            }
        }
    }
    for (Chain* C : binderChains) {
        for (Residue* R : C->getResidues()) {
            // only store alpha carbon atoms
            Atom* CA = R->findAtom("CA",true);
            binder.push_back(CA);
        }
    }

    AtomPointerVector aVec(target);
    ProximitySearch PS(aVec,cutoff/2);
    int count = 0;
    for (Atom* CA : binder) {
        Residue* binderR = CA->getResidue();
        set<Residue*> nearbyResidues;
        vector<int> nearbyAtoms = PS.getPointsWithin(CA->getCoor(),0,cutoff);
        for (int atomIdx : nearbyAtoms) {
            Residue* R = target[atomIdx]->getResidue();
            nearbyResidues.insert(R);
        }
        for (Residue* R : nearbyResidues) {
            interfaceResidues.push_back(pair<Residue*,Residue*>(R,binderR));
            count++;
        }
    }
    // cout << "Found " << count << " interface contacts" << endl;
}

residueFrame* interfaceScorer::defineQuery(Residue* reference, Residue* mobile) {
    residueFrame* refFrame = complex->getResidueFrame(reference->getResidueIndex());
    residueFrame* mobileFrame = complex->getResidueFrame(mobile->getResidueIndex());
    return mobileFrame->frameRelativeToOther(*refFrame);
}