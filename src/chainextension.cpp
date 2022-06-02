#include "chainextension.h"

/* --- --- --- --- --- extensionTopology --- --- --- --- --- */

Structure& extensionTopology::getTerminalSegment(extensionDirection terminus) {
    if (terminus == extensionDirection::CTERM) {
        mappedSegment& ctermSegment = allSegments.back();
        return ctermSegment.segment;
    } else {
        mappedSegment& ntermSegment = allSegments.front();
        return ntermSegment.segment;
    }
}

void extensionTopology::addSegment(const Structure& segment, extensionDirection extDir, int numResExt, bool fixed) {
    // Make a copy of the structure object that will be managed by this class
    int ntermPos;
    if (extDir == extensionDirection::NTERM) {
        mappedSegment& ntermSegment = allSegments.front();
        ntermPos = ntermSegment.ntermPos - numResExt;
        // mappedSegment mSegment(ntermPos,segmentCopy,fixed);
        allSegments.emplace_front(ntermPos,segment,fixed);
    } else {
        mappedSegment& ctermSegment = allSegments.back();
        ntermPos = ctermSegment.ctermPos() + numResExt - (segment.residueSize() - 1);
        // mappedSegment mSegment(ntermPos,segmentCopy,fixed);
        allSegments.emplace_back(ntermPos,segment,fixed);
    }
    fusedStructureCurrent = false;
};

bool extensionTopology::removeSegment(extensionDirection terminus) {
    if (terminus == extensionDirection::CTERM) {
        mappedSegment& ctermSegment = allSegments.back();
        if (ctermSegment.fixed) return false; // cannot remove fixed segment
        allSegments.pop_back();
        return true;
    } else {
        mappedSegment& ntermSegment = allSegments.front();
        if (ntermSegment.fixed) return false; // cannot remove fixed segment
        allSegments.pop_front();
        return true;
    }
    fusedStructureCurrent = false;
};

fusionTopology extensionTopology::buildTopology() {
    // Find total length and starting point
    int minPos = allSegments.front().ntermPos;
    int maxPos = allSegments.back().ctermPos();
    int totalLength = maxPos - minPos + 1;
    // cout << "minPos: " << minPos << " maxPos: " << maxPos << " totalLength: " << totalLength << endl;

    fusionTopology fTopology(totalLength);
    for (mappedSegment& mSeg : allSegments) {
        vector<int> segResIdx(mSeg.numRes);
        int absoluteNtermPos = mSeg.ntermPos-minPos;
        // cout << "absoluteNtermPos: " << absoluteNtermPos << endl;
        iota(begin(segResIdx),end(segResIdx),absoluteNtermPos);
        // for (int i : segResIdx) cout << i << " ";
        // cout << endl;
        fTopology.addFragment(mSeg.segment,segResIdx);
        // if (mSeg.fixed) fTopology.addFixedPositions(segResIdx);
    }
    return fTopology;
};

void extensionTopology::writeTopologyToPDB(string namePrefix) {
    string name = namePrefix + ".topology.pdb";
    fstream pdb_out;
    MstUtils::openFile(pdb_out,name,fstream::out);
    for (mappedSegment& mSeg : allSegments) {
        pdb_out << "HEADER    " << namePrefix << "_" << mSeg.ntermPos << "_" << mSeg.numRes << endl;
        mSeg.segment.writePDB(pdb_out);
    }
}

Structure& extensionTopology::getFusedStructure() {
    if (fusedStructureCurrent) return fused;

    /* Borrowed the following fuser code from TERMify
    fuses first without internal coordinate constraints
    after aligning all the fragments, then focuses on satisfying the constraints
    */
    vector<mstreal> ics = fParams.getIntCoorFCs();
    fParams.setIntCoorFCs({0, 0, 0});
    fParams.clearStartingStructure();
    fused = Fuser::fuse(buildTopology(),fParams);
    fParams.setIntCoorFCs(ics);
    fParams.setStartingStructure(fused);
    fused = Fuser::fuse(buildTopology(),fParams);
    fusedStructureCurrent = true;
    return fused;
};

void extensionTopology::writeFusedStructureToPDB(string name) {
    if (!fusedStructureCurrent) getFusedStructure();
    fused.writePDB(name+".pdb");
}

void extensionTopology::writeFusedStructureToPDB(fstream& out, string name, int state) {
    out << "MODEL state" << state << "_" << name << endl;
    if (!fusedStructureCurrent) getFusedStructure();
    fused.writePDB(out);
    out << "ENDMDL" << endl;
}

int extensionTopology::getResidueSize() {
    if (!fusedStructureCurrent) getFusedStructure();
    return fused.residueSize();
}

/* --- --- --- --- --- sampleSegmentExtensionsFASST --- --- --- --- --- */

sampleSegmentExtensionsFASST::sampleSegmentExtensionsFASST(string fasstDBPath) {
    cout << "Loading DB..." << endl;
    F.readDatabase(fasstDBPath,1);
    cout << "Done loading DB." << endl;
};

int sampleSegmentExtensionsFASST::searchForMatchesToFragment(Structure& querySeg, int numResExt, extensionDirection extDir) {
    // Find matches to the query
    F.setQuery(querySeg);
    allMatches = F.search();

    // Select the matches that contain sufficient flanking residues
    viableMatches.clear();
    for (set<fasstSolution>::iterator it = allMatches.begin(); it != allMatches.end(); it++) attemptToExtend(*it,numResExt,extDir);
    indexSampler = fisherYatesShuffle(0,viableMatches.size()-1);
    return viableMatches.size();
};

Structure* sampleSegmentExtensionsFASST::sampleExtensionFragment() {
    int currentMatchIndex = indexSampler.sample();
    if (currentMatchIndex < 0) MstUtils::error("Must search for matches before sampling extension fragments");
    if (currentMatchIndex > viableMatches.size() - 1) MstUtils::error("No match corresponding to the index");
    Structure* segment = getSegment(viableMatches[currentMatchIndex]);
    return segment;
};

void sampleSegmentExtensionsFASST::attemptToExtend(const fasstSolution& sol, int numResExt, extensionDirection extDir) {
    int targetIdx, alignment, segLength;
    targetIdx = sol.getTargetIndex();
    alignment = sol[0];
    segLength = sol.segLength(0);
    if ((targetIdx < 0)||(targetIdx >= F.numTargets())) MstUtils::error("Target with index: "+MstUtils::toString(targetIdx)+" does not exist");
    
    Structure* target = F.getTarget(targetIdx);
    Residue* alignRes = &target->getResidue(alignment);
    Chain* C = alignRes->getChain();
    int alignResIdx = alignRes->getResidueIndex();
    int alignResIdxInChain = alignRes->getResidueIndexInChain();

    // Get the residue range (inclusive) in the chain
    int NtermResIdxInChain, CtermResIdxInChain;
    if (extDir == extensionDirection::NTERM) {
        NtermResIdxInChain = alignResIdxInChain - numResExt;
        CtermResIdxInChain = alignResIdxInChain + (segLength - 1);
    } else {
        NtermResIdxInChain = alignResIdxInChain;
        CtermResIdxInChain = alignResIdxInChain + (segLength - 1) + numResExt;
    }
    if ((NtermResIdxInChain < 0)||(CtermResIdxInChain >= C->residueSize())) return; // cannot be extended

    // if sufficient flanking residues to extend the fragment, store the info for later
    int newAlignment = NtermResIdxInChain + (alignResIdx - alignResIdxInChain);
    int totalLength = segLength + numResExt;
    matchWithExtension extSol(sol,newAlignment,totalLength);
    viableMatches.push_back(extSol);
};

Structure* sampleSegmentExtensionsFASST::getSegment(matchWithExtension match) {
    int targetIdx = match.sol.getTargetIndex();
    if ((targetIdx < 0)||(targetIdx >= F.numTargets())) MstUtils::error("Target with index: "+MstUtils::toString(targetIdx)+" does not exist");
    Structure* target = F.getTarget(targetIdx);
    vector<Residue*> selectedTargetResidues;
    for (int i = match.alignment; i < match.alignment + match.totalLength; i++) selectedTargetResidues.push_back(&target->getResidue(i));
    Structure* segment = new Structure(selectedTargetResidues);
    Transform tf = match.sol.getTransform();
    tf.apply(segment);
    return segment;
}

/* --- --- --- --- --- binderChainExtensionFASST --- --- --- --- --- */

void binderChainExtensionFASST::setFixedStructuralContext(Structure& context, const binderScorerParams& sParams) {
    target = augmentedStructure(context);
    scorer = new binderScorer(sParams,target);
    scorer->defineTargetBindingSiteResiduesByrSASA();
    targetBB = Structure(MiscTools::getBackboneAtoms(target));
    clashCheck.setStructure(targetBB);
};

void binderChainExtensionFASST::extendChainRandom(int totalExtensionLength) {
    if (!readyToExtend()) MstUtils::error("Object not prepared before extension");

    // Before extending, score and save initial structure
    int round = 0;
    binder = augmentedStructure(extTopology.getFusedStructure());
    openTrajectoryFile("binderExtension.traj.pdb");
    saveStructureToTrajectoryFile(binder,MstUtils::toString(round));

    round++;
    int currentExtensionLength = 0, initialStructureLength = binder.residueSize(), maxAttempts = 10;
    while (currentExtensionLength < totalExtensionLength) {
        extensionDirection currentTerminus;
        if (terminus != extensionDirection::EITHER) currentTerminus = terminus;
        else currentTerminus = (MstUtils::randInt(2)) ? extensionDirection::NTERM : extensionDirection::CTERM;

        Structure terminalSegment = getTerminalSegment(currentTerminus);
        sampler.searchForMatchesToFragment(terminalSegment,extensionSegmentLength,terminus);
        cout << "Found " << sampler.numRemainingSegments() << " viable matches" << endl;

        // find a valid extension fragment 
        int attemptCount = 0;
        Structure* extensionSegment = nullptr;
        Structure binderWithExtension;
        while ((attemptCount < maxAttempts)&(sampler.numRemainingSegments() > 0)) {
            cout << "Round: " << round << " attemptCount: " << attemptCount << endl;
            extensionSegment = sampler.sampleExtensionFragment();
            extensionSegment->writePDB("extensionSegment_"+MstUtils::toString(round)+"_"+MstUtils::toString(attemptCount)+".pdb");
            
            // add to topology and fuse
            extTopology.addSegment(*extensionSegment,currentTerminus,extensionSegmentLength);
            delete extensionSegment;
            binderWithExtension = extTopology.getFusedStructure();

            // check for clashes
            if (!extensionSegmentClash(&binderWithExtension)) {
                // clash between the extension segment and the protein
                break;
            }

            saveStructureToTrajectoryFile(binderWithExtension,MstUtils::toString(round)+"_"+MstUtils::toString(attemptCount));
            extTopology.removeSegment(currentTerminus);
            cout << "extension fragment clashed, sampling another... (" << sampler.numRemainingSegments() << " remaining)" << endl;
            attemptCount++;
        }
        if (attemptCount >= maxAttempts) {
            cout << "Could not find an extension fragment that did not clash. Removing terminal segment and finding new extension..." << endl;
            extTopology.removeSegment(currentTerminus);
            binder = augmentedStructure(extTopology.getFusedStructure());
            continue;
        }
        binder = binderWithExtension;

        saveStructureToTrajectoryFile(binder,MstUtils::toString(round)+"_"+MstUtils::toString(attemptCount));
        extTopology.writeTopologyToPDB("round_"+MstUtils::toString(round));

        currentExtensionLength = binder.residueSize() - initialStructureLength;
        cout << "currentExtensionLength: " << currentExtensionLength << endl;
        round++;
    }
    traj_out.close();
};

void binderChainExtensionFASST::extendChainGreedy(int totalExtensionLength, int NOptionsPerStep) {
    if (!readyToExtend()) MstUtils::error("Object not prepared before extension");

    // Before extending, score and save initial structure
    int round = 0, step = 1;
    binder = augmentedStructure(extTopology.getFusedStructure());
    openTrajectoryFile("binderExtension.traj.pdb");
    saveStructureToTrajectoryFile(binder,MstUtils::toString(round));

    fstream data_out;
    MstUtils::openFile(data_out,"binderExtensionScores.csv",fstream::out);
    data_out << "step,name,numRes,numPotConts,numNonDesignable,score" << endl;

    scorer->setBinder(&binder);
    scorer->defineInterfaceByPotentialContacts();
    mstreal score = scorer->scoreInterface();
    cout << "score: " << score << endl;
    cout << "nondesignable contacts: " << scorer->countNonDesignableContacts() << endl;
    data_out << step << ",";
    data_out << "0_1," << binder.residueSize() << "," << scorer->countDesignableContacts() << ",";
    data_out << scorer->countNonDesignableContacts() << "," << score << endl;

    round++;
    int currentExtensionLength = 0, initialStructureLength = binder.residueSize(), maxAttempts = 10;
    while (currentExtensionLength < totalExtensionLength) {
        extensionDirection currentTerminus;
        if (terminus != extensionDirection::EITHER) currentTerminus = terminus;
        else currentTerminus = (MstUtils::randInt(2)) ? extensionDirection::NTERM : extensionDirection::CTERM;

        Structure terminalSegment = getTerminalSegment(currentTerminus);
        sampler.searchForMatchesToFragment(terminalSegment,extensionSegmentLength,terminus);
        cout << "Found " << sampler.numRemainingSegments() << " viable matches" << endl;

        // Find N valid extensions and select the one with the best score
        int potentialExtensionCount = 1;
        Structure* proposedExtensionSegment = nullptr;
        Structure* bestExtensionSegment = nullptr;
        mstreal bestExtensionSegmentTotalScore = 0.0;
        augmentedStructure binderWithExtension;
        while ((potentialExtensionCount <= NOptionsPerStep)&(sampler.numRemainingSegments() > 0)) {
            proposedExtensionSegment = sampler.sampleExtensionFragment();
            // extensionSegment->writePDB("extensionSegment_"+MstUtils::toString(round)+"_"+MstUtils::toString(attemptCount)+".pdb");
            
            // add to topology and fuse
            extTopology.addSegment(*proposedExtensionSegment,currentTerminus,extensionSegmentLength);
            binderWithExtension = augmentedStructure(extTopology.getFusedStructure());

            if (extensionSegmentClash(&binderWithExtension)) {
                // if clashing, find another extension
                extTopology.removeSegment(currentTerminus);
                delete proposedExtensionSegment;
                cout << "extension fragment clashed, sampling another... (" << sampler.numRemainingSegments() << " remaining)" << endl;
                continue;
            }
            // otherwise, score and consider for extension
            scorer->setBinder(&binderWithExtension);
            scorer->defineInterfaceByPotentialContacts();
            score = scorer->scoreInterface();
            cout << "score: " << score << endl;
            cout << "nondesignable contacts: " << scorer->countNonDesignableContacts() << endl;
            data_out << step << ",";
            data_out << MstUtils::toString(round)+"_"+MstUtils::toString(potentialExtensionCount) << ",";
            data_out << binderWithExtension.residueSize() << ",";
            data_out << scorer->countDesignableContacts() << ",";
            data_out << scorer->countNonDesignableContacts() << "," << score << endl;

            // set as new best
            if (score < bestExtensionSegmentTotalScore) {
                bestExtensionSegmentTotalScore = score;
                delete bestExtensionSegment;
                bestExtensionSegment = proposedExtensionSegment;
            } else {
                delete proposedExtensionSegment;
            }

            saveStructureToTrajectoryFile(binderWithExtension,MstUtils::toString(round)+"_"+MstUtils::toString(potentialExtensionCount));
            extTopology.removeSegment(currentTerminus);

            potentialExtensionCount++;
            step++;
        }

        // Select the best of the valid extensions and repeat
        cout << "End of round " << round << " selecting extension with score " << bestExtensionSegmentTotalScore << endl;
        extTopology.addSegment(*bestExtensionSegment,currentTerminus,extensionSegmentLength);
        delete bestExtensionSegment;
        binder = augmentedStructure(extTopology.getFusedStructure());
        extTopology.writeTopologyToPDB("round_"+MstUtils::toString(round));

        currentExtensionLength = binder.residueSize() - initialStructureLength;
        cout << "currentExtensionLength: " << currentExtensionLength << endl;
        round++;
    }
    binder.writePDB("extendedBinder.fin.pdb");
    data_out.close();
    traj_out.close();
};

bool binderChainExtensionFASST::readyToExtend() {
    if (anchorSegment.residueSize() == 0) return false;
    if (target.residueSize() == 0) return false;
    return true;
};

Structure binderChainExtensionFASST::getTerminalSegment(extensionDirection terminus, bool allowShorter) {
    int startResIdx, endResIdx; //[startResIdx, endResIdx)
    int currentOverlapSegmentLength = overlapSegmentLength;
    if (binder.residueSize() < overlapSegmentLength) {
        if (allowShorter) currentOverlapSegmentLength = binder.residueSize();
        else MstUtils::error("Starting segment shorter than the desired overlap length","binderChainExtensionFASST::getTerminalSegment");
    }
    if (terminus == extensionDirection::NTERM) {
        startResIdx = 0;
        endResIdx = currentOverlapSegmentLength;
    } else {
        startResIdx = binder.residueSize() - currentOverlapSegmentLength;
        endResIdx = binder.residueSize();
    }
    vector<Residue*> allRes = binder.getResidues();
    return Structure(vector<Residue*>(allRes.begin()+startResIdx,allRes.begin()+endResIdx));
};

bool binderChainExtensionFASST::extensionSegmentClash(Structure* extensionSegment) {
    // check for clashes between the binder and target
    if (clashCheck.checkForClashesToStructure(extensionSegment->getResidues())) return true;
    // check for clashes within the binder
    if (clashCheck.checkForClashesWithinQueryStructure(extensionSegment->getResidues())) return true;
    return false;
}

void binderChainExtensionFASST::openTrajectoryFile(string pathToFile) {
    MstUtils::openFile(traj_out,pathToFile,fstream::out);
};

void binderChainExtensionFASST::saveStructureToTrajectoryFile(const Structure& S, string stateName) {
    traj_out << "MODEL " << stateName << endl;
    S.writePDB(traj_out);
    traj_out << "ENDMDL" << endl;
};


/* --- --- --- --- --- binderChainExtension --- --- --- --- --- */

binderChainExtension::binderChainExtension(Structure& context, string overlapGraphPath, const binderScorerParams& sParams) : overlaps(overlapGraphPath) {
    // set target
    target = augmentedStructure(context);
    scorer = new binderScorer(sParams,target);
    scorer->defineTargetBindingSiteResiduesByrSASA();
    targetBB = Structure(MiscTools::getBackboneAtoms(target));
    clashCheck.setStructure(targetBB);
};

void binderChainExtension::setBinderAnchor(Structure& anchor) {
    // set binder seed
    anchorSegment = anchor;
    anchorName = MstSys::splitPath(anchor.getName(),1);
}

void binderChainExtension::extendChain(int lengthToExtend, extensionDirection terminus, extensionMode mode, int beamWidth) {
    if (anchorSegment.residueSize() == 0) MstUtils::error("Must set the binder anchor before extension","binderChainExtension::extendChain");

    // Initialize the set of binder extensions with the single anchor
    vector<scoredBinder*> selectedBinders;
    vector<scoredBinder*> selectedBinderExtensions;

    openTrajectoryFile(anchorName+".traj.pdb");
    state_count = 1;

    scoredBinder* binder = new scoredBinder(anchorSegment,"anchor",state_count);
    binder->topology.writeFusedStructureToPDB(traj_out,binder->name,state_count);
    state_count++;

    mstreal score = scoreBinder(*binder,0);
    cout << "Starting score: " << score << endl;
    selectedBinders.push_back(binder);

    int anchorLength = anchorSegment.residueSize();
    int extensionLength = binder->topology.getResidueSize() - anchorLength;
    extensionDirection terminusToExtendThisRound;
    int round = 1;
    while (extensionLength < lengthToExtend) {
        cout << "Round: " << round << endl;

        // Select the terminus to extend
        if (terminus != extensionDirection::EITHER) terminusToExtendThisRound = terminus;
        else terminusToExtendThisRound = (MstUtils::randInt(2)) ? extensionDirection::NTERM : extensionDirection::CTERM;

        // Generate all viable extensions of the current topologies
        for (scoredBinder* currentBinder : selectedBinders) {
            Structure* terminalSegment = &currentBinder->topology.getTerminalSegment(terminusToExtendThisRound);
            timer.start();
            vector<pair<Structure*,int>> extensionSegments = overlaps.getExtensionSegments(terminalSegment,terminusToExtendThisRound);
            timer.stop();
            cout << "Found " << extensionSegments.size() << " potential extensions to " << currentBinder->getName() << " in " << timer.getDuration() << " s" << endl;
            int extensionNumber = 1;
            for (pair<Structure*,int> extensionSegment: extensionSegments) {
                cout << extensionSegment.first->getName() << endl;
                if (extensionSegmentClash(extensionSegment.first)) {
                    // clash between the binder chain and the protein
                    cout << "Skipping due to clash between binder extension segment and the protein..." << endl;
                    continue;
                }
                scoredBinder* newBinder = new scoredBinder(*currentBinder);
                newBinder->name+="_"+MstUtils::toString(round)+"-"+MstUtils::toString(extensionNumber);
                newBinder->topology.addSegment(*extensionSegment.first,terminusToExtendThisRound,overlaps.getOverlapLength(),false);
                
                if (binderIntrachainClash(newBinder)) {
                    // clash between the binder chain and itself
                    cout << "Skipping due to intra-chain clash" << endl;
                    delete newBinder;
                    continue;
                }

                cout << "Added " << extensionSegment.first->getName() << " with " << extensionSegment.first->getResidues().size() << " residues" << endl;
                score = scoreBinder(*newBinder,extensionSegment.second);
                cout << "Score: " << score << " name: " << newBinder->name << endl;
                selectedBinderExtensions.push_back(newBinder);
                newBinder->topology.writeFusedStructureToPDB(traj_out,newBinder->name,state_count);
                state_count++;
                extensionNumber++;
            }
        }
        if (selectedBinderExtensions.empty()) {
            cout << "Could not find any viable extensions, terminating..." << endl;
            break;
        }

        // Sort by score
        cout << "Sorting " << selectedBinderExtensions.size() << " potential extensions by score" << endl;
        timer.start();
        // for (scoredBinder* binder : selectedBinderExtensions) cout << binder->score << " ";
        // cout << endl;
        if (mode != extensionMode::RANDOM) sort(selectedBinderExtensions.begin(),selectedBinderExtensions.end(),[](scoredBinder* a, scoredBinder* b) { return a->score < b->score; });
        // for (scoredBinder* binder : selectedBinderExtensions) cout << binder->score << " ";
        // cout << endl;
        timer.stop();
        cout << "Took " << timer.getDuration() << " to sort the selected binder extensions" << endl;

        // Select the top candidate(s) according to the mode
        for (scoredBinder* binder : selectedBinders) delete binder;
        selectedBinders.clear();
        if (mode == extensionMode::RANDOM) {
            int randomSel = MstUtils::randInt(selectedBinderExtensions.size());
            scoredBinder* sel = selectedBinderExtensions[randomSel];
            cout << "Selecting " << sel->name << " state (" << state_count << ") with score = " << sel->score << endl;
            for (int i = 0; i < selectedBinderExtensions.size(); i++) {
                if (i == randomSel) selectedBinders.push_back(selectedBinderExtensions[i]);
                else delete selectedBinderExtensions[i];
            }
        } else if (mode == extensionMode::GREEDY) {
            scoredBinder* sel = selectedBinderExtensions.front();
            cout << "Selecting " << sel->name << " with score = " << sel->score << endl;
            selectedBinders.push_back(sel);
            for (int i = 1; i < selectedBinderExtensions.size(); i++) delete selectedBinderExtensions[i];
        } else if (mode == extensionMode::BEAM) {
            for (int i = 0; i < selectedBinderExtensions.size(); i++) {
                if (i < beamWidth) {
                    scoredBinder* sel = selectedBinderExtensions[i];
                    cout << "Selecting " << sel->name << " with score = " << sel->score << endl;
                    selectedBinders.push_back(sel);
                }
                else delete selectedBinderExtensions[i];
            }
        }
        selectedBinderExtensions.clear();
        cout << "Selected " << selectedBinders.size() << " binders for further extension" << endl;

        extensionLength = selectedBinders.front()->topology.getResidueSize() - anchorLength;
        round++;
    }
    // Write out selected binders to PDB files
    for (scoredBinder* binder : selectedBinders) {
        binder->topology.writeFusedStructureToPDB(binder->name+".fin");
        delete binder;
    }
    closeTrajectoryFile();
};

void binderChainExtension::coverChainWithExtensionSegments(Chain* chainToCover, extensionDirection terminus, int beamWidth) {
    if (terminus == extensionDirection::EITHER) MstUtils::error("Must either select NTERM or CTERM","binderChainExtension::coverChainWithExtensionSegments");

    // Define the anchor segment according to the provided peptide
    vector<Atom*> chainBBAtoms = MiscTools::getBackboneAtoms(chainToCover);
    Structure StoCover(chainBBAtoms);
    vector<Residue*> Sresidues = StoCover.getResidues();
    vector<Residue*> anchorRes;
    if (terminus == extensionDirection::CTERM) anchorRes = vector<Residue*>(Sresidues.begin(),Sresidues.begin()+overlaps.getOverlapLength());
    else anchorRes = vector<Residue*>(Sresidues.end()-overlaps.getOverlapLength(),Sresidues.end());
    anchorSegment = Structure(anchorRes);
    anchorName = chainToCover->getParent()->getName();
    
    openTrajectoryFile(anchorName+".traj.pdb");
    state_count = 1;

    vector<scoredBinder*> binderExtensions;
    vector<scoredBinder*> selectedBinders;
    scoredBinder* binder = new scoredBinder(anchorSegment,"anchor",state_count);
    selectedBinders.push_back(binder);
    binder->topology.writeFusedStructureToPDB(traj_out,binder->name,state_count);
    state_count++;

    // mstreal score = scoreBinder(*binder,0);
    // cout << "Starting score: " << score << endl;

    int anchorLength = anchorSegment.residueSize();
    int extensionLength = binder->topology.getResidueSize() - anchorLength;
    int lengthToExtend = chainToCover->residueSize() - anchorLength;
    int round = 1;
    vector<Residue*> chainSection;
    vector<Atom*> chainSectionBBAtoms;
    RMSDCalculator calc;
    while (extensionLength < lengthToExtend - overlaps.getOverlapLength() + 1) {
        // Get the section of the chain that is to be covered
        if (terminus == extensionDirection::CTERM) chainSection = vector<Residue*>(Sresidues.begin()+(round*overlaps.getOverlapLength()),Sresidues.begin()+((round+1)*overlaps.getOverlapLength()));
        else chainSection = vector<Residue*>(Sresidues.end()-((round+1)*overlaps.getOverlapLength()),Sresidues.end()-(round*overlaps.getOverlapLength()));
        chainSectionBBAtoms = MiscTools::getBackboneAtoms(chainSection);

        cout << "Round: " << round << endl;

        // Generate all viable extensions of the current topologies
        for (scoredBinder* currentBinder : selectedBinders) {
            // Generate all viable extensions of the current topologies
            Structure* terminalSegment = &currentBinder->topology.getTerminalSegment(terminus);
            timer.start();
            vector<pair<Structure*,int>> extensionSegments = overlaps.getExtensionSegments(terminalSegment,terminus);
            timer.stop();
            cout << "Found " << extensionSegments.size() << " potential extensions to " << currentBinder->getName() << " in " << timer.getDuration() << " s" << endl;
            int extensionNumber = 1;
            for (pair<Structure*,int> extensionSegment: extensionSegments) {
                cout << extensionSegment.first->getName() << endl;
                if (extensionSegmentClash(extensionSegment.first)) {
                    // clash between the binder chain and the protein
                    cout << "Skipping due to clash between binder extension segment and the protein..." << endl;
                    continue;
                }
                scoredBinder* newBinder = new scoredBinder(*currentBinder);
                newBinder->name+="_"+MstUtils::toString(round)+"-"+MstUtils::toString(extensionNumber);
                newBinder->topology.addSegment(*extensionSegment.first,terminus,overlaps.getOverlapLength(),false);
                
                if (binderIntrachainClash(newBinder)) {
                    // clash between the binder chain and itself
                    cout << "Skipping due to intra-chain clash" << endl;
                    delete newBinder;
                    continue;
                }

                cout << "Added " << extensionSegment.first->getName() << " with " << extensionSegment.first->getResidues().size() << " residues" << endl;
                // score = scoreBinder(*newBinder,extensionSegment.second);
                // cout << "Score: " << score << " name: " << newBinder->name << endl;

                newBinder->topology.writeFusedStructureToPDB(traj_out,newBinder->name,state_count);

                // Replace the score value with the RMSD to the chainToCover
                vector<Residue*> extensionResidues = newBinder->topology.getFusedStructure().getResidues();
                vector<Residue*> extensionSection;
                if (terminus == extensionDirection::CTERM) extensionSection = vector<Residue*>(extensionResidues.end()-overlaps.getOverlapLength(),extensionResidues.end());
                else extensionSection = vector<Residue*>(extensionResidues.begin(),extensionResidues.begin()+overlaps.getOverlapLength());
                vector<Atom*> extensionSectionBBAtoms = MiscTools::getBackboneAtoms(extensionSection);
                newBinder->score = calc.rmsd(chainSectionBBAtoms,extensionSectionBBAtoms);
                cout << "RMSD: " << newBinder->score << endl;

                binderExtensions.push_back(newBinder);
                state_count++;
                extensionNumber++;
            }
        }
        if (binderExtensions.empty()) {
            cout << "Could not find any viable extensions, terminating..." << endl;
            break;
        }

        // Sort by RMSD to the chain that is to be covered
        cout << "Sorting " << binderExtensions.size() << " potential extensions by RMSD" << endl;
        timer.start();
        // for (scoredBinder* binder : selectedBinderExtensions) cout << binder->score << " ";
        // cout << endl;s
        sort(binderExtensions.begin(),binderExtensions.end(),[](scoredBinder* a, scoredBinder* b) { return a->score < b->score; });
        // for (scoredBinder* binder : selectedBinderExtensions) cout << binder->score << " ";
        // cout << endl;
        timer.stop();
        cout << "Took " << timer.getDuration() << " to sort the binder extensions" << endl;

        // Select the top candidates
        for (scoredBinder* binder : selectedBinders) delete binder;
        selectedBinders.clear();
        for (int i = 0; i < binderExtensions.size(); i++) {
            if (i < beamWidth) {
                scoredBinder* sel = binderExtensions[i];
                cout << "Selecting " << sel->name << " with RMSD = " << sel->score << endl;
                selectedBinders.push_back(sel);
            }
            else delete binderExtensions[i];
        }
        binderExtensions.clear();
        
        extensionLength = selectedBinders.front()->topology.getResidueSize() - anchorLength;
        round++;
    }
    // Write out selected binder as a PDB file
    for (scoredBinder* binder : selectedBinders) {
        binder->topology.writeFusedStructureToPDB(binder->name+".fin");
        delete binder;
    }
    closeTrajectoryFile();
}

bool binderChainExtension::extensionSegmentClash(Structure* extensionSegment) {
    // check for clashes between the binder and target
    if (clashCheck.checkForClashesToStructure(extensionSegment->getResidues())) return true;
    // check for clashes within the binder
    if (clashCheck.checkForClashesWithinQueryStructure(extensionSegment->getResidues())) return true;
    return false;
}

bool binderChainExtension::binderIntrachainClash(scoredBinder* binder) {
    return clashCheck.checkForClashesWithinQueryStructure(binder->getFusedStructure().getResidues());
}

mstreal binderChainExtension::scoreBinder(scoredBinder& topology, int segmentDesignability) {
    augmentedStructure binder = topology.getFusedStructure();
    scorer->setBinder(&binder);
    scorer->defineInterfaceByPotentialContacts();
    mstreal score = scorer->scoreInterface();
    score += (-log(segmentDesignability+1)*overlaps.getOverlapLength());
    topology.setScore(score);
    return score;
}

void binderChainExtension::openTrajectoryFile(string pathToFile) {
    MstUtils::openFile(traj_out,pathToFile,fstream::out);
};