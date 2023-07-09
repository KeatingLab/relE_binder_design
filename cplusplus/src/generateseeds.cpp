#include "generateseeds.h"

/* --- --- --- --- --- seedGenerator --- --- --- --- --- */

seedGenerator::seedGenerator(const Structure& _target, string fasstDBPath, seedGenParams _params) : target(_target), params(_params) {
    targetName = MstSys::splitPath(target.getName(),1);

    defineBindingSiteAsSurfaceRes();

    for (Residue* R : target.getResidues()) {
        if (find(bindingSite.begin(),bindingSite.end(),R) == bindingSite.end()) {
            vector<Atom*> resAtoms = R->getAtoms();
            targetBackboneAtoms.insert(targetBackboneAtoms.end(),resAtoms.begin(),resAtoms.end());
        } else if (RotamerLibrary::hasFullBackbone(R)) {
            vector<Atom*> resBBAtoms = RotamerLibrary::getBackbone(R);
            targetBackboneAtoms.insert(targetBackboneAtoms.end(),resBBAtoms.begin(),resBBAtoms.end());
        }
    }
    targetAtomsPS = ProximitySearch(AtomPointerVector(targetBackboneAtoms),1.5);

    if (params.contactData != "") {
        potConts.load2DProbabilityDensities(params.contactData);
        potConts.setTargetResidues(target.getResidues());
    }

    cout << "Load FASSTDB (" << fasstDBPath << ") without sidechains" << endl;
    timer.start();
    F.readDatabase(fasstDBPath,1);
    timer.stop();
    cout << "Done took " << timer.getDuration() << " seconds to load DB" << endl;

    // set 
}

void seedGenerator::defineBindingSiteAsSurfaceRes(mstreal rSASAthresh) {
    bindingSite.empty();

    sasaCalculator calc(target);
    calc.computeSASA();
    map<Residue*,mstreal> resSASA = calc.getResidueSASA(true);
    for (auto it : resSASA) {
        if (it.second > rSASAthresh) bindingSite.push_back(it.first);
    }
    cout << "Found " << bindingSite.size() << " surface residues with a relative SASA threshold of " << rSASAthresh << endl;
}

string seedGenerator::generateSeeds() {
    if (bindingSite.empty()) MstUtils::error("Cannot generate seeds, no binding site residues defined","seedGenerator::generateSeeds()");

    if (!allFragmentsData.empty()) for (auto frag : allFragmentsData) delete frag.fragment;
    allFragmentsData.clear();

    // open seed binary / data files
    string seedBinPath = targetName+".seeds.bin";
    seedBinaryFile seedBin(seedBinPath,false);
    fstream fragOut, seedOut;
    MstUtils::openFile(fragOut,targetName+"_fragments.csv",fstream::out);
    fragOut << "fragName,cenResIdx,chainID,resNum,numMatches,numSeeds" << endl;
    MstUtils::openFile(seedOut,targetName+"_seeds.csv",fstream::out);
    seedOut << "seedName,numRes,target_idx,target_name,match_res_chain_id,match_res_num,seed_res_chain_id,seed_res_num" << endl;
    multiPDBFile* pdbOut = nullptr;
    if (params.writeToPDB) {
        pdbOut = new multiPDBFile(targetName+".seeds.pdb",false);
        pdbOut->addStructure(target,targetName+"_target");
    }

    // define a fragment around each binding site residue
    getBindingSiteFragments();

    // search each fragment against the DB
    timer.start();
    findStructuralMatches();

    // generate seed(s) from each match, if possible
    generateSeedsFromMatches(seedBin,fragOut,seedOut,pdbOut);

    fragOut.close(); seedOut.close();
    return seedBinPath;
}

void seedGenerator::writeBindingSiteFragments(string path) {
    if (allFragmentsData.empty()) MstUtils::error("No fragments to write to file","seedGenerator::writeBindingSiteFragments()");

    fstream out;
    MstUtils::openFile(out,path,fstream::out);

    for (int i = 0; i < allFragmentsData.size(); i++) {
        fragmentData& fD = allFragmentsData[i];
        out << "HEADER    " << fD.getName() << endl;
        fD.fragment->writePDB(out);
    }
}

void seedGenerator::getBindingSiteFragments() {
    if (bindingSite.empty()) MstUtils::error("Cannot generate seeds without first defining a binding site");
    for (Residue* R : bindingSite) {
        fragmentData fD;
        fD.parent = this;
        fD.fragment = new Structure();
        fD.cenResIdxInFrag = TERMUtils::selectTERM({R},*fD.fragment,params.targetFlankRes)[0];
        fD.fragCenResInParent = R;
        allFragmentsData.push_back(fD);
    }
}

void seedGenerator::findStructuralMatches() {
    for (fragmentData& fD : allFragmentsData) {
        // set query and params
        F.setQuery(*fD.fragment);
        F.setRMSDCutoff(params.RMSDCutoff);
        F.setMaxNumMatches(params.maxNumMatches);
        F.options().unsetSequenceConstraints();

        if (params.seqConst) {
            cout << "Set a sequence constraint before searching the fragment..." << endl;
            fasstSeqConstSimple seqConst(1);
            seqConst.addConstraint(0,fD.cenResIdxInFrag,{fD.fragCenResInParent->getName()});
            F.options().setSequenceConstraints(seqConst);
        }

        // search
        cout << "Searching fragment: " << fD.getName() << endl;
        timer.start();
        fD.matches = F.search();
        timer.stop();
        cout << "Took " << timer.getDuration() << " seconds to find " << fD.matches.size() << " matches" << endl;
        // vector<Sequence> matchSeq = F.getMatchSequences(fD.matches);
        // for (Sequence s : matchSeq) cout << s.toString() << endl;
    }
}

void seedGenerator::generateSeedsFromMatches(seedBinaryFile& seedBin, fstream& fragOut, fstream& seedOut, multiPDBFile* pdbOut) {
    for (fragmentData& fD: allFragmentsData) {
        cout << "Generate seeds from fragment: " << fD.getName() << endl;
        int seedCount = 0;
        for (fasstSolution sol : fD.matches) {
            // get residues contacting the central residue of the match
            vector<int> matchIndices = F.getMatchResidueIndices(sol);
            Residue* cenMatchRes = &F.getTarget(sol.getTargetIndex())->getResidue(matchIndices[fD.cenResIdxInFrag]);
            if (!F.isResiduePairBoolPropertyDefined(sol.getTargetIndex(),"vdw")) MstUtils::error("vdw property not defined for this structure in DB","seedGenerator::generateSeedsFromMatches");
            set<int> contactingRes = F.getResiduePairBoolProperty(sol.getTargetIndex(),"vdw",matchIndices[fD.cenResIdxInFrag]);

            if (contactingRes.empty()) continue;

            // for each contact, make a seed
            Structure* target = F.getTarget(sol.getTargetIndex());
            map<Structure*,Residue*> seed2contactR;
            vector<Structure*> seeds;
            for (int contactingResIdx : contactingRes) {
                // use chain boundaries to restrict which flanking residues are included in the seed
                Residue* contactR = &target->getResidue(contactingResIdx);
                Chain* matchChain = contactR->getChain();
                int resIdxInChain = matchChain->getResidueIndex(contactR);
                int startInChain = max(0,resIdxInChain - params.seedFlankRes);
                int endInChain = min(matchChain->residueSize() - 1,resIdxInChain + params.seedFlankRes);
                
                // convert index in chain to index in structure
                int start = startInChain + (contactingResIdx - resIdxInChain);
                int end = endInChain + (contactingResIdx - resIdxInChain);

                vector<Residue*> seedRes;
                Structure* seed = new Structure(); seeds.push_back(seed);
                // fragmentname_targetidx_matchresidx_seedresidx
                // string seedName = fD.getName() + "_" + MstUtils::toString(sol.getTargetIndex()) + "_" 
                // + MstUtils::toString(matchIndices[fD.cenResIdxInFrag]) + "_" + MstUtils::toString(contactingResIdx);
                string seedName = getSeedName(sol.getTargetIndex(),start,end-start+1);
                seed->setName(seedName);

                Chain* seedChain = seed->appendChain("0",false);
                int resNum = 1;
                for (int i = start; i <= end; i++) {
                    Residue* targetR = &target->getResidue(i);
                    if (!RotamerLibrary::hasFullBackbone(targetR)) MstUtils::error("Residue "+MstUtils::toString(targetR->getResidueIndex())+" from target PDB "+MstUtils::toString(sol.getTargetIndex())+" is missing a backbone atom","seedGenerator::generateSeedsFromMatches");
                    Residue* Rcopy = new Residue(target->getResidue(i));
                    Rcopy->setIcode(' ');
                    Rcopy->setNum(resNum);
                    if (!SeqToolsExtension::AAinSet(Rcopy->getName())) Rcopy->setName("ALA");
                    seedChain->appendResidue(Rcopy);
                    resNum++;
                }

                seed2contactR[seed] = contactR;
            }

            // transform each seed and see if it clashes with the target backbone
            Transform tr = sol.getTransform();
            vector<Structure*> transformedSeeds;
            for (Structure* seed : seeds) {
                tr.apply(seed);

                if (seedTargetClash(seed)) {
                    delete seed;
                    continue;
                }
                transformedSeeds.push_back(seed);
            }

            // add the anchor chain
            // TO DO (if needed for fusing)

            // store
            for (Structure* seed : transformedSeeds) {
                if ((params.contactData != "")&&(params.minContsPerRes > 0)) {
                    potConts.setBinderResidues(seed->getResidues());
                    potConts.getContacts();
                    int nConts = potConts.getNumContacts();
                    mstreal contsPerRes = mstreal(nConts)/mstreal(seed->residueSize());
                    if (contsPerRes < params.minContsPerRes) {
                        delete seed;
                        continue;
                    } 
                }
                seedBin.appendStructure(seed);

                seedOut << seed->getName() << "," << seed->residueSize() << ",";
                seedOut << sol.getTargetIndex() << "," << F.getTargetName(sol.getTargetIndex()) << ",";
                seedOut << cenMatchRes->getChainID() << "," << cenMatchRes->getNum() << ",";
                seedOut << seed2contactR[seed]->getChainID() << "," << seed2contactR[seed]->getNum();
                seedOut << endl;

                if (pdbOut != nullptr) pdbOut->addStructure(*seed);

                delete seed;
                seedCount++;
            }
        }
        fragOut << fD.getName() << "," << fD.fragCenResInParent->getChainID() << ",";
        fragOut << fD.fragCenResInParent->getNum() << "," << fD.matches.size() << ",";
        fragOut << seedCount << endl;
    }
}

bool seedGenerator::seedTargetClash(Structure* seed) {
    int checkDistance = checker.maxSumRadii();
    for (Residue* R : seed->getChainByID("0")->getResidues()) {
        for (Atom* A : R->getAtoms()) {
            vector<int> nearbyAtomIdx = targetAtomsPS.getPointsWithin(A->getCoor(),0,checkDistance);
            // now check if any of these nearby atoms are bona-fide clashes
            for (int atomIdx : nearbyAtomIdx) {
                if (checker.clash(A,*targetBackboneAtoms[atomIdx])) return true;
            }
        }
    }
    return false;
}

