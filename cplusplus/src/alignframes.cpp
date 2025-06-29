#include "alignframes.h"

/* --- --- --- --- --- proteinFrameDB --- --- --- --- --- */

const augmentedStructure& augmentedStructureDB::getTarget(int i) {
    if ((i < 0) || (i >= targets.size()))
        MstUtils::error("Provided value " + MstUtils::toString(i) + " is out of range: (0," + MstUtils::toString(targets.size() - 1), "proteinFrameDB::getTarget");
    return *targets[i];
}

Residue *augmentedStructureDB::getResidue(int target_i, int res_i) {
    if ((target_i < 0) || (target_i >= targets.size())) MstUtils::error("Provided value " + MstUtils::toString(target_i) + " is out of range of targets: (0," + MstUtils::toString(targets.size() - 1), "proteinFrameDB::getResidue");
    if ((res_i < 0) || (res_i >= targets[target_i]->residueSize())) MstUtils::error("Provided value " + MstUtils::toString(res_i) + " is out of range of residues: (0," + MstUtils::toString(targets[target_i]->residueSize() - 1), "proteinFrameDB::getResidue");
    return &targets[target_i]->getResidue(res_i);
}

residueFrame *augmentedStructureDB::getResidueFrame(int target_i, int res_i) {
    if ((target_i < 0) || (target_i >= targets.size())) MstUtils::error("Provided value " + MstUtils::toString(target_i) + " is out of range of targets: (0," + MstUtils::toString(targets.size() - 1), "proteinFrameDB::getResidueFrame");
    if ((res_i < 0) || (res_i >= targets[target_i]->residueSize())) MstUtils::error("Provided value " + MstUtils::toString(res_i) + " is out of range of residues: (0," + MstUtils::toString(targets[target_i]->residueSize() - 1), "proteinFrameDB::getResidueFrame");
    return targets[target_i]->getResidueFrame(res_i);
}

void augmentedStructureDB::addTarget(augmentedStructure *S)
{
    targets.push_back(S);
}

void augmentedStructureDB::addTarget(const augmentedStructure &S)
{
    augmentedStructure *newTarget = new augmentedStructure(S);
    targets.push_back(newTarget);
}

void augmentedStructureDB::readDBFile(string dbPath, int debug_N)
{
    cout << "Reading the structure DB..." << endl;
    MstTimer timer; timer.start();
    fstream ifs;
    MstUtils::openFile(ifs, dbPath, fstream::in | fstream::binary, "alignInteractingFrames::readDBFile");
    
    char sect;
    string name;
    int ver = 0;
    int ti = numTargets();
    MstUtils::readBin(ifs, sect);
    if (sect == 'V')
    {
        MstUtils::readBin(ifs, ver);
        MstUtils::readBin(ifs, sect);
    }
    if (sect != 'S')
        MstUtils::error("first section must be a structure one, while reading database file " + dbPath, "alignInteractingFrames::readDBFile");
    while (ifs.peek() != EOF)
    {
        Structure targetStruct;
        streampos loc = ifs.tellg();
        targetStruct.readData(ifs);
        augmentedStructure *targetAugS = new augmentedStructure(targetStruct);
        // targetSource.push_back(targetInfo(dbFile, targetFileType::BINDATABASE, loc, memSave));
        int L = targetStruct.residueSize();
        addTarget(targetAugS);
        while (ifs.peek() != EOF) {
            MstUtils::readBin(ifs, sect);
            if (sect == 'B') {
                MstUtils::readBin(ifs, name);
                if ((name != "vdw")&&(name != "vdwBB")) MstUtils::error("Section with name: "+name+" not recognized","proteinFrameDB::readDBFile");
                map<int, set<int>> &vals = (name == "vdw") ? vdwContacts[ti] : vdwContactsBB[ti];
                int ri, rj, N, n;
                MstUtils::readBin(ifs, N);
                for (int i = 0; i < N; i++) {
                    MstUtils::readBin(ifs, ri);
                    MstUtils::readBin(ifs, n);
                    for (int j = 0; j < n; j++) {
                        MstUtils::readBin(ifs, rj);
                        vals[ri].insert(rj);
                    }
                }
            }
            else if (sect == 'S')
            {
                break;
            }
            else
            {
                MstUtils::error("unknown section type: " + MstUtils::toString(sect) + ", while reading database file " + dbPath, "alignInteractingFrames::readDBFile");
            }
        }
        ti++;
        if ((debug_N > 0)&&(ti >= debug_N)) break;
    }
    ifs.close();
    timer.stop();
    cout << "Loaded database with " << targets.size() << " structures in " << timer.getDuration() << " s" << endl;
}

void augmentedStructureDB::writeDBFile(string dbPath)
{
    fstream ofs;
    MstUtils::openFile(ofs, dbPath, fstream::out | fstream::binary, "alignInteractingFrames::writeDBFile");
    MstUtils::writeBin(ofs, 'V');
    MstUtils::writeBin(ofs, (int)version); // format version
    for (int ti = 0; ti < targets.size(); ti++)
    {
        if (targets[ti] == NULL) MstUtils::error("Cannot write a database in which full structures are not populated", "alignInteractingFrames::writeDBFile");
        MstUtils::writeBin(ofs, 'S'); // marks the start of a structure section
        targets[ti]->writeData(ofs);

        if (vdwContacts.find(ti) != vdwContacts.end())
        {
            map<int, set<int>> &vals = vdwContacts[ti];
            MstUtils::writeBin(ofs, 'B'); // marks the start of a residue pair interaction property section
            MstUtils::writeBin(ofs, (string)"vdw");
            MstUtils::writeBin(ofs, (int)vals.size());
            for (auto i = vals.begin(); i != vals.end(); ++i)
            {
                MstUtils::writeBin(ofs, (int)i->first);
                MstUtils::writeBin(ofs, (int)(i->second).size());
                for (int j : i->second) 
                {
                    MstUtils::writeBin(ofs, (int)j);
                }
            }
        }
    }
    ofs.close();
}

const set<int>& augmentedStructureDB::getContacts(int target_i, int res_i, bool BBint) {
    if ((target_i < 0) || (target_i >= targets.size())) MstUtils::error("Provided value " + MstUtils::toString(target_i) + " is out of range of targets: (0," + MstUtils::toString(targets.size() - 1), "proteinFrameDB::getContacts");
    if ((res_i < 0) || (res_i >= targets[target_i]->residueSize())) MstUtils::error("Provided value " + MstUtils::toString(res_i) + " is out of range of residues: (0," + MstUtils::toString(targets[target_i]->residueSize() - 1), "proteinFrameDB::getContacts");
    return vdwContacts[target_i][res_i];
}

/* --- --- --- --- --- alignInteractingFrames --- --- --- --- --- */

void alignInteractingFrames::setAA(string resName) {
    aa = SeqTools::aaToIdx(resName);
    aaName = SeqTools::idxToTriple(aa);
}

void alignInteractingFrames::findMobileFrames() {
    if (subsample_flanking > 0) {
        cout << "Flanking residue subsample rate set to " << subsample_flanking << endl;
        MstUtils::seedRandEngine(42);
    }
    // if (aa == SeqTools::unknownIdx())
    //     MstUtils::error("Must set amino acid before finding interacting residues", "alignInteractingFrames::findInteractingRes");
    // allInteractingFrames.clear();
    allInteractingRes.clear();
    // interactionData[aa] = map<int,int>();
    // map<int,int>& interactionDataForAA = interactionData[aa];
    // for (int i = 0; i <= 30; i++) interactionDataForAA[i] = 0;
    int cis_count = 0;
    int cic_count = 0;
    for (int target_idx = 0; target_idx < db.numTargets(); target_idx++) {
        const augmentedStructure& aS = db.getTarget(target_idx);
        PC->setTargetResidues(aS.getResidues());
        PC->setBinderResidues(aS.getResidues());
        if (verbose) cout << "target " << target_idx << " has " << aS.residueSize() << " residues and " << aS.chainSize() << " chains" << endl;
        for (int res_i = 0; res_i < aS.residueSize(); res_i++) {
            Residue* Ri = &aS.getResidue(res_i);
            res_t res_i_aa = SeqTools::aaToIdx(Ri->getName());
            if ((aa != SeqTools::unknownIdx())&(res_i_aa != aa)) {
                continue;
            }

            // Find residues that contact res_i
            set<int> contacting_res;
            if (PC == nullptr) {
                contacting_res = db.getContacts(target_idx,res_i);
            } else {
                const set<Residue*> contacting_res_pointers = PC->getContactsWithResidue(Ri);
                for (Residue* R : contacting_res_pointers) contacting_res.insert(R->getResidueIndex());
            }
            // interactionDataForAA[contacting_res.size()]++;
            for (int res_j : contacting_res) {
                Residue* Rj = &aS.getResidue(res_j);
                int distance_in_chain = -1; // interaction via contact (long range)
                allInteractingRes.emplace_back(Ri,Rj,target_idx,distance_in_chain);
                cis_count++;
                // mobileFrame* frame = new mobileFrame(Rj, target_idx, res_i, res_j, res_i_aa);
                // residueFrame *rFi = db.getResidueFrame(frame->getTarget(),frame->getResI());
                // Transform rFi_to_ref = TransformFactory::switchFrames(refFrame, *rFi);
                // rFi_to_ref.apply(*frame);
                // allInteractingFrames.push_back(frame);
                if (verbose) cout << target_idx << "," << res_i << "," << res_j << endl;
            }

            // Find residues that are proximal to res_i in the chain
            if (flanking_res <= 0) continue;
            int chain_start_res_idx = res_i - Ri->getResidueIndexInChain();
            int chain_stop_res_idx = chain_start_res_idx + Ri->getChain()->residueSize() - 1;
            if (verbose) cout << "chain start: " << chain_start_res_idx << " chain stop: " << chain_stop_res_idx << endl;
            for (int res_j = max(chain_start_res_idx,res_i - flanking_res); res_j <= min(chain_stop_res_idx,res_i + flanking_res); res_j++) {
                if (res_i == res_j) continue;
                // randomly skip some flanking residues in order to keep the total number from exploding
                if (subsample_flanking >= 0) if (MstUtils::randUnit() > subsample_flanking) continue;
                Residue* Rj = &aS.getResidue(res_j);
                int distance_in_chain = abs(res_j - res_i);
                allInteractingRes.emplace_back(Ri,Rj,target_idx,distance_in_chain);
                cic_count++;
                // mobileFrame* frame = new mobileFrame(Rj, target_idx, res_i, res_j, res_i_aa);
                // residueFrame *rFi = db.getResidueFrame(frame->getTarget(),frame->getResI());
                // Transform rFi_to_ref = TransformFactory::switchFrames(refFrame, *rFi);
                // rFi_to_ref.apply(*frame);
                // allInteractingFrames.push_back(frame);
                if (verbose) cout << target_idx << "," << res_i << "," << res_j << endl;
            }
        }
    }
    cout << "In total, found " << cis_count << " close-in-space contact res pairs and " << cic_count << " close-in-chain contact res pairs" << endl;
}

// Structure *alignInteractingFrames::getAlignedInteractingRes(int i) {
//     if ((i < 0) || (i >= allInteractingFrames.size())) MstUtils::error("Provided value " + MstUtils::toString(i) + " is out of range: (0," + MstUtils::toString(allInteractingFrames.size() - 1), "alignInteractingFrames::getAlignedInteractingRes");
//     mobileFrame* frame = allInteractingFrames[i];

//     return getAlignedInteractingRes(frame);
// };

// Structure *alignInteractingFrames::getAlignedInteractingRes(mobileFrame* frame) {
//     // get the residues from the target structure
//     Residue* rI = db.getResidue(frame->getTarget(),frame->getResI());
//     Residue* rJ = db.getResidue(frame->getTarget(),frame->getResJ());

//     Structure *interactingResiduePair = constructStructureFromResiduePair(rI, rJ, frame);

//     return changeFrameToRef(interactingResiduePair,frame);
// }

// void alignInteractingFrames::writeAlignedInteractingResToPDB(string pdbPath, mstreal subsampleRate)
// {
//     // open file
//     fstream pdb_out;
//     MstUtils::openFile(pdb_out, pdbPath, fstream::out, "alignInteractingFrames::writeAlignedInteractingResToPDB");

//     // transform interacting residues
//     for (int i = 0; i < getNumInteracting(); i++)
//     {
//         if (MstUtils::randUnit() > subsampleRate) continue;
//         Structure *interactingResiduePair = getAlignedInteractingRes(i);

//         // write to file
//         pdb_out << "HEADER    " << interactingResiduePair->getName() << endl;
//         interactingResiduePair->writePDB(pdb_out);

//         delete interactingResiduePair;
//     }
//     pdb_out.close();
// }

// void alignInteractingFrames::writeMobileFramesToBin(frameDB* frameBin) {
//     for (mobileFrame* mF : allInteractingFrames) {
//         frameBin->appendFrame(mF);
//     }
// }

void alignInteractingFrames::writeResiduePairsToBin(resPairDB* rPBin) {
    for (int i = 0; i < allInteractingRes.size(); i++) {
        rPBin->appendResPair(&allInteractingRes[i]);
    }
}

void alignInteractingFrames::writeInteractionData(string pathPrefix) {
    fstream info_out;
    string path = pathPrefix + "_interactionData.tsv";
    MstUtils::openFile(info_out,path,fstream::out,"alignInteractingFrames::writeInteractionData");
    info_out << "aminoAcidType\tnumLongRangeContacts\tcount" << endl;

    for (auto it_aa : interactionData) {
        for (auto it_cont : it_aa.second) {
            info_out << SeqTools::idxToTriple(it_aa.first) << "\t";
            info_out << it_cont.first << "\t";
            info_out << it_cont.second;
            info_out << endl;
        }
    }
    info_out.close();
}

bool alignInteractingFrames::isQueryHomologousToMatchInDB(Residue* query, Residue* reference, mobileFrame* mFrame) {
    Residue* qRi = query;
    Residue* qRj = reference;
    Residue* mRi = db.getResidue(mFrame->getTarget(),mFrame->getResI());
    Residue* mRj = db.getResidue(mFrame->getTarget(),mFrame->getResJ());

    pair<int,int> Ri_identity = getSequenceIdentity(qRi,mRi);
    pair<int,int> Rj_identity = getSequenceIdentity(qRj,mRj);

    mstreal seqID = mstreal(Ri_identity.second + Rj_identity.second) / mstreal(Ri_identity.first + Rj_identity.first);
    // cout << "seqID: " << seqID << " over " << (Ri_identity.first + Rj_identity.first) << " residues" << endl;
    return seqID >= homologyThreshold;
}

Structure *alignInteractingFrames::changeFrameToRef(Structure* interactingResiduePair, mobileFrame* frame) {
    // get the transformation that aligns residue i to the reference frame
    residueFrame *rFi = db.getResidueFrame(frame->getTarget(),frame->getResI());
    Transform rFi_to_ref = TransformFactory::switchFrames(refFrame, *rFi);
    // Transform rFi_to_ref = TransformFactory::switchFrames(*rFi,rFrame);

    // apply the transformation to the selected residues
    rFi_to_ref.apply(interactingResiduePair);
    return interactingResiduePair;
}

Structure *alignInteractingFrames::constructStructureFromResiduePair(Residue *Ri, Residue *Rj, mobileFrame* frame)
{
    Structure *interactingResiduePair = new Structure;

    // add each residue
    Chain *chainA = new Chain("A", "");
    chainA->insertResidueCopy(Ri);
    bool allow_rename = false;
    interactingResiduePair->appendChain(chainA, allow_rename);

    Chain *chainB = new Chain("B", "");
    chainB->insertResidueCopy(Rj);
    interactingResiduePair->appendChain(chainB, allow_rename);

    string name = frame->getName();
    interactingResiduePair->setName(name);
    return interactingResiduePair;
}

pair<int,int> alignInteractingFrames::getSequenceIdentity(Residue* R1, Residue* R2) {
    Chain* R1Chain = R1->getChain();
    Chain* R2Chain = R2->getChain();

    if ((R1Chain == NULL)||(R1Chain == NULL)) MstUtils::error("Somehow residues in the database have no parent chains","alignInteractingFrames::getSequenceIdentity");

    int R1Index = R1->getResidueIndexInChain();
    int R2Index = R2->getResidueIndexInChain();

    // Adjust residue indices so that R1 and R2 have idx 0, so that ranges can be compared
    int R1_min = max(R1Index-windowSize,0) - R1Index;
    int R1_max = min(R1Index+windowSize,R1Chain->residueSize()-1) - R1Index;
    int R2_min = max(R2Index-windowSize,0) - R2Index;
    int R2_max = min(R2Index+windowSize,R2Chain->residueSize()-1) - R2Index;

    int windowMin = max(R1_min,R2_min);
    int windowMax = min(R1_max,R2_max);

    // Extract the sequence windows from the complete chain sequence
    Sequence R1ChainSeq(*R1Chain);
    Sequence R2ChainSeq(*R2Chain);

    Sequence R1WindowSeq = R1ChainSeq.extractRange(R1Index+windowMin,R1Index+windowMax);
    Sequence R2WindowSeq = R2ChainSeq.extractRange(R2Index+windowMin,R2Index+windowMax);

    if (R1WindowSeq.length() != R1WindowSeq.length()) MstUtils::error("Sequences must be the same length","alignInteractingFrames::getSequenceIdentity");
    return pair<int,int>(R1WindowSeq.length(),SeqTools::sequenceIdentity(R1WindowSeq,R2WindowSeq));
}

/* --- --- --- --- --- mobileFrame --- --- --- --- --- */

void mobileFrame::writeData(ostream& ofs) {
    // write each point
    writeCartesianPointToBin(getO(),ofs);
    writeCartesianPointToBin(getX(),ofs);
    writeCartesianPointToBin(getY(),ofs);
    writeCartesianPointToBin(getZ(),ofs);

    // write DB location
    MstUtils::writeBin(ofs,target);
    MstUtils::writeBin(ofs,res_i);
    MstUtils::writeBin(ofs,res_j);

    // write res types
    MstUtils::writeBin(ofs,res_i_aa);
    MstUtils::writeBin(ofs,res_j_aa);
}

void mobileFrame::readData(istream& ifs) {
    CartesianPoint O, X, Y, Z;
    readCartesianPointFromBin(O,ifs);
    readCartesianPointFromBin(X,ifs);
    readCartesianPointFromBin(Y,ifs);
    readCartesianPointFromBin(Z,ifs);

    MstUtils::readBin(ifs, target);

    MstUtils::readBin(ifs, res_i);
    MstUtils::readBin(ifs, res_j);

    MstUtils::readBin(ifs, res_i_aa);
    MstUtils::readBin(ifs, res_j_aa);

    constructFrame(O,X,Y,Z);
}

void mobileFrame::writeCartesianPointToBin(const CartesianPoint& p, ostream& ofs) {
    MstUtils::writeBin(ofs,p.getX());
    MstUtils::writeBin(ofs,p.getY());
    MstUtils::writeBin(ofs,p.getZ());
}

void mobileFrame::readCartesianPointFromBin(CartesianPoint& p, istream& ifs) {
    p.resize(3);
    mstreal& x = p[0];
    mstreal& y = p[1];
    mstreal& z = p[2];
    MstUtils::readBin(ifs,x);
    MstUtils::readBin(ifs,y);
    MstUtils::readBin(ifs,z);
}

/* --- --- --- --- --- frameDB --- --- --- --- --- */

bool frameDB::hasNext() {
    if (!readMode) MstUtils::error("hasNext not supported in write mode","frameDB::hasNext");
    return fs.peek() != EOF;
}

void frameDB::skip() {
    if (!readMode) MstUtils::error("skip not supported in write mode","frameDB::skip");
    mobileFrame* frame = readNextFileSection();
    delete frame;
}

void frameDB::reset() {
    if (!readMode) MstUtils::error("reset not supported in write mode","frameDB::reset");
    fs.clear();
    fs.seekg(0,fs.beg);
}

mobileFrame* frameDB::next() {
    if (!readMode) MstUtils::error("next not supported in write mode","frameDB::next");
    return readNextFileSection();
}

vector<mobileFrame*> frameDB::loadAllFrames() {
    vector<mobileFrame*> loaded; 
    reset();
    while (hasNext()) {
        mobileFrame* mF = next();
        loaded.push_back(mF);
    }
    return loaded;
}

void frameDB::appendFrame(mobileFrame* frame) {
    if (readMode) MstUtils::error("appendStructure not supported in write mode","frameDB::appendStructure");
    if (!frameAdded) {
        frameAdded = true;
        MstUtils::writeBin(fs,version);
    }
    MstUtils::writeBin(fs,'F'); //start new frame section
    frame->writeData(fs);
    MstUtils::writeBin(fs,'E'); //the E is not necessary, but gives me flexibility in case I ever want to store more data in the file
}

void frameDB::openFileStream() {
    if (readMode) MstUtils::openFile(fs, frameDBPath, fstream::in | fstream::binary, "frameDB::openFileStream");
    else MstUtils::openFile(fs, frameDBPath, fstream::out | fstream::binary, "frameDB::openFileStream");
}

mobileFrame* frameDB::readNextFileSection() {
    //if beginning of file, advance past the version
    if (fs.tellg() == 0) {
        int version; MstUtils::readBin(fs, version);
    }
    char sect;
    MstUtils::readBin(fs, sect);
    if (sect != 'F') MstUtils::error("The first section should be a residueFrame. frame database: " + frameDBPath, "frameDB::readNextFileSection");
    mobileFrame* frame = new mobileFrame;
    frame->readData(fs);
    MstUtils::readBin(fs, sect);
    if (sect != 'E') MstUtils::error("The data should terminate with 'E'. frame database: " + frameDBPath, "frameDB::readNextFileSection");

    return frame;
}

/* --- --- --- --- --- resPair --- --- --- --- --- */

resPair::resPair(Residue* Ri, Residue* Rj, int _target, int _distance_in_chain) : target(_target), distance_in_chain(_distance_in_chain) {
    res_i = Ri->getResidueIndex();
    res_j = Rj->getResidueIndex();
    res_i_aa = SeqTools::aaToIdx(Ri->getName());
    res_j_aa = SeqTools::aaToIdx(Rj->getName());
    vector<Atom*> RiBBAtoms = RotamerLibrary::getBackbone(Ri);
    vector<Atom*> RjBBAtoms = RotamerLibrary::getBackbone(Rj);
    resPairAtoms.insert(resPairAtoms.end(),RiBBAtoms.begin(),RiBBAtoms.end());
    resPairAtoms.insert(resPairAtoms.end(),RjBBAtoms.begin(),RjBBAtoms.end());
    // computeDistancesFromBBAtoms;
    computeInternalRepresentation();
}

// void resPair::computeDistancesFromBBAtoms() {
//     mstreal NDistance = resPairAtoms[0]->getCoor().distance(resPairAtoms[4]->getCoor());
//     mstreal CaDistance = resPairAtoms[1]->getCoor().distance(resPairAtoms[5]->getCoor());
//     mstreal CDistance = resPairAtoms[2]->getCoor().distance(resPairAtoms[6]->getCoor());
//     bbAtomDistances = CartesianPoint(NDistance,CaDistance,CDistance);
// }

CartesianPoint resPair::getAllbbAtomDistances() {
    // Compute all inter-residue distances between the 4 backbone atoms of the residue pair
    mstreal N_N   = resPairAtoms[0]->getCoor().distance(resPairAtoms[4]->getCoor());
    mstreal N_Ca  = resPairAtoms[0]->getCoor().distance(resPairAtoms[5]->getCoor());
    mstreal N_C   = resPairAtoms[0]->getCoor().distance(resPairAtoms[6]->getCoor());
    mstreal N_O   = resPairAtoms[0]->getCoor().distance(resPairAtoms[7]->getCoor());
    mstreal Ca_N  = resPairAtoms[1]->getCoor().distance(resPairAtoms[4]->getCoor());
    mstreal Ca_Ca = resPairAtoms[1]->getCoor().distance(resPairAtoms[5]->getCoor());
    mstreal Ca_C  = resPairAtoms[1]->getCoor().distance(resPairAtoms[6]->getCoor());
    mstreal Ca_O  = resPairAtoms[1]->getCoor().distance(resPairAtoms[7]->getCoor());
    mstreal C_N   = resPairAtoms[2]->getCoor().distance(resPairAtoms[4]->getCoor());
    mstreal C_Ca  = resPairAtoms[2]->getCoor().distance(resPairAtoms[5]->getCoor());
    mstreal C_C   = resPairAtoms[2]->getCoor().distance(resPairAtoms[6]->getCoor());
    mstreal C_O   = resPairAtoms[2]->getCoor().distance(resPairAtoms[7]->getCoor());
    mstreal O_N   = resPairAtoms[3]->getCoor().distance(resPairAtoms[4]->getCoor());
    mstreal O_Ca  = resPairAtoms[3]->getCoor().distance(resPairAtoms[5]->getCoor());
    mstreal O_C   = resPairAtoms[3]->getCoor().distance(resPairAtoms[6]->getCoor());
    mstreal O_O   = resPairAtoms[3]->getCoor().distance(resPairAtoms[7]->getCoor());
    return CartesianPoint({N_N,N_Ca,N_C,N_O,
                          Ca_N,Ca_Ca,Ca_C,Ca_O,
                          C_N,C_Ca,C_C,C_O,
                          O_N,O_Ca,O_C,O_O});
}

void resPair::computeInternalRepresentation() {
    CaDistance = resPairAtoms[1]->getCoor().distance(resPairAtoms[5]->getCoor());
    residueFrame Ri(resPairAtoms[0],resPairAtoms[1],resPairAtoms[2]);
    residueFrame Rj(resPairAtoms[4],resPairAtoms[5],resPairAtoms[6]);

    // Get all three basis vectors for each residue frame
    CartesianPoint Ri_x = Ri.getX();
    CartesianPoint Rj_x = Rj.getX();
    CartesianPoint Ri_y = Ri.getY();
    CartesianPoint Rj_y = Rj.getY();
    CartesianPoint Ri_z = Ri.getZ();
    CartesianPoint Rj_z = Rj.getZ();

    // // Compute the cosine angle between each pair of basis vectors
    // mstreal xCosSim = Ri_x.cosineAngle(Rj_x);
    // mstreal yCosSim = Ri_y.cosineAngle(Rj_y);
    // mstreal zCosSim = Ri_z.cosineAngle(Rj_z);

    // mstreal xAngle, yAngle, zAngle;
    // if (xCosSim == 1.0) xAngle = 0.0;
    // else xAngle = (180/M_PI)*acos(xCosSim);
    // if (yCosSim == 1.0) yAngle = 0.0;
    // else yAngle = (180/M_PI)*acos(yCosSim);
    // if (zCosSim == 1.0) zAngle = 0.0;
    // else zAngle = (180/M_PI)*acos(zCosSim);
    // resFrameBasisVectorAngles = CartesianPoint(xAngle,yAngle,zAngle);
    // if (std::isnan(xAngle) || std::isnan(yAngle) || std::isnan(zAngle)) {
    //     cout << "Ca distance: " << getCaDistance() << endl;
    //     cout << "Ri_x: " << Ri_x << endl;
    //     cout << "Rj_x: " << Rj_x << endl;
    //     cout << "Ri_y: " << Ri_y << endl;
    //     cout << "Rj_y: " << Rj_y << endl;
    //     cout << "Ri_z: " << Ri_z << endl;
    //     cout << "Rj_z: " << Rj_z << endl;
    //     cout << "xCosSimilarity: " << Ri_x.cosineAngle(Rj_x) << endl;
    //     cout << "xAngle radians: " << acos(Ri_x.cosineAngle(Rj_x)) << endl;
    //     cout << "yCosSimilarity: " << Ri_y.cosineAngle(Rj_y) << endl;
    //     cout << "yAngle radians: " << acos(Ri_y.cosineAngle(Rj_y)) << endl;
    //     cout << "zCosSimilarity: " << Ri_z.cosineAngle(Rj_z) << endl;
    //     cout << "zAngle radians: " << acos(Ri_z.cosineAngle(Rj_z)) << endl;
    // }

    // Get the distances between backbone atoms
    mstreal NDistance = resPairAtoms[0]->getCoor().distance(resPairAtoms[4]->getCoor());
    mstreal CaDistance = resPairAtoms[1]->getCoor().distance(resPairAtoms[5]->getCoor());
    mstreal CDistance = resPairAtoms[2]->getCoor().distance(resPairAtoms[6]->getCoor());
    bbAtomDistances = CartesianPoint(NDistance,CaDistance,CDistance);
}

void resPair::writeData(ostream& ofs) {
    // write each atom
    for (int i = 0; i < resPairAtoms.size(); i++) writeAtomToBin(resPairAtoms[i],ofs);

    // write DB location
    MstUtils::writeBin(ofs,target);
    MstUtils::writeBin(ofs,res_i);
    MstUtils::writeBin(ofs,res_j);

    MstUtils::writeBin(ofs,distance_in_chain);

    // write res types
    MstUtils::writeBin(ofs,res_i_aa);
    MstUtils::writeBin(ofs,res_j_aa);
}

void resPair::readData(istream& ifs) {
    for (int i = 0; i < 8; i++) {
        Atom* A = readAtomFromBin(ifs);
        resPairAtoms.push_back(A);
    }
    ownsAtoms = true;

    MstUtils::readBin(ifs, target);
    MstUtils::readBin(ifs, res_i);
    MstUtils::readBin(ifs, res_j);

    MstUtils::readBin(ifs, distance_in_chain);

    MstUtils::readBin(ifs, res_i_aa);
    MstUtils::readBin(ifs, res_j_aa);

    computeInternalRepresentation();
}

void resPair::writeAtomToBin(Atom* A, ostream& ofs) {
    // write x,y,x coordinate (discard the rest of the information)
    MstUtils::writeBin(ofs,A->getX());
    MstUtils::writeBin(ofs,A->getY());
    MstUtils::writeBin(ofs,A->getZ());
}

Atom* resPair::readAtomFromBin(istream& ifs) {
    Atom* A = new Atom();
    mstreal x, y, z;
    MstUtils::readBin(ifs,x);
    MstUtils::readBin(ifs,y);
    MstUtils::readBin(ifs,z);
    A->setX(x);
    A->setY(y);
    A->setZ(z);
    A->stripInfo();
    return A;
}

/* --- --- --- --- --- resPairDB --- --- --- --- --- */

bool resPairDB::hasNext() {
    if (!readMode) MstUtils::error("hasNext not supported in write mode","resPairDB::hasNext");
    return fs.peek() != EOF;
}

void resPairDB::skip() {
    if (!readMode) MstUtils::error("skip not supported in write mode","resPairDB::skip");
    resPair* rP = readNextFileSection();
    delete rP;
}

void resPairDB::reset() {
    if (!readMode) MstUtils::error("reset not supported in write mode","resPairDB::reset");
    fs.clear();
    fs.seekg(0,fs.beg);
}

resPair* resPairDB::next() {
    if (!readMode) MstUtils::error("next not supported in write mode","resPairDB::next");
    return readNextFileSection();
}

vector<resPair*> resPairDB::loadAllResPairs() {
    vector<resPair*> loaded; 
    reset();
    while (hasNext()) {
        resPair* mF = next();
        loaded.push_back(mF);
    }
    return loaded;
}

void resPairDB::appendResPair(resPair* rP) {
    if (readMode) MstUtils::error("appendResPair not supported in read mode","resPairDB::appendStructure");
    if (!resPairAdded) {
        resPairAdded = true;
        MstUtils::writeBin(fs,version);
    }
    MstUtils::writeBin(fs,'F'); //start new frame section
    rP->writeData(fs);
    MstUtils::writeBin(fs,'E'); //the E is not necessary, but gives me flexibility in case I ever want to store more data in the file
}

void resPairDB::openFileStream() {
    if (readMode) MstUtils::openFile(fs, resPairDBPath, fstream::in | fstream::binary, "resPairDB::openFileStream");
    else MstUtils::openFile(fs, resPairDBPath, fstream::out | fstream::binary, "resPairDB::openFileStream");
}

resPair* resPairDB::readNextFileSection() {
    //if beginning of file, advance past the version
    if (fs.tellg() == 0) {
        int version; MstUtils::readBin(fs, version);
    }
    char sect;
    MstUtils::readBin(fs, sect);
    if (sect != 'F') MstUtils::error("The first section should be a residueFrame. frame database: " + resPairDBPath, "resPairDB::readNextFileSection");
    resPair* rP = new resPair;
    rP->readData(fs);
    MstUtils::readBin(fs, sect);
    if (sect != 'E') MstUtils::error("The data should terminate with 'E'. frame database: " + resPairDBPath, "resPairDB::readNextFileSection");

    return rP;
}