#include "alignframes.h"

/* --- --- --- --- --- proteinFrameDB --- --- --- --- --- */

const augmentedStructure& proteinFrameDB::getTarget(int i) {
    if ((i < 0) || (i >= targets.size()))
        MstUtils::error("Provided value " + MstUtils::toString(i) + " is out of range: (0," + MstUtils::toString(targets.size() - 1), "proteinFrameDB::getTarget");
    return *targets[i];
}

Residue *proteinFrameDB::getResidue(int target_i, int res_i) {
    if ((target_i < 0) || (target_i >= targets.size())) MstUtils::error("Provided value " + MstUtils::toString(target_i) + " is out of range of targets: (0," + MstUtils::toString(targets.size() - 1), "proteinFrameDB::getResidue");
    if ((res_i < 0) || (res_i >= targets[target_i]->residueSize())) MstUtils::error("Provided value " + MstUtils::toString(res_i) + " is out of range of residues: (0," + MstUtils::toString(targets[target_i]->residueSize() - 1), "proteinFrameDB::getResidue");
    return &targets[target_i]->getResidue(res_i);
}

residueFrame *proteinFrameDB::getResidueFrame(int target_i, int res_i) {
    if ((target_i < 0) || (target_i >= targets.size())) MstUtils::error("Provided value " + MstUtils::toString(target_i) + " is out of range of targets: (0," + MstUtils::toString(targets.size() - 1), "proteinFrameDB::getResidueFrame");
    if ((res_i < 0) || (res_i >= targets[target_i]->residueSize())) MstUtils::error("Provided value " + MstUtils::toString(res_i) + " is out of range of residues: (0," + MstUtils::toString(targets[target_i]->residueSize() - 1), "proteinFrameDB::getResidueFrame");
    return targets[target_i]->getResidueFrame(res_i);
}

void proteinFrameDB::addTarget(augmentedStructure *S)
{
    targets.push_back(S);
}

void proteinFrameDB::addTarget(const augmentedStructure &S)
{
    augmentedStructure *newTarget = new augmentedStructure(S);
    targets.push_back(newTarget);
}

void proteinFrameDB::readDBFile(string dbPath)
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
                MstUtils::error("unknown section type" + MstUtils::toString(sect) + ", while reading database file " + dbPath, "alignInteractingFrames::readDBFile");
            }
        }
        ti++;
    }
    ifs.close();
    timer.stop();
    cout << "Loaded database with " << targets.size() << " structures in " << timer.getDuration() << " s" << endl;
}

void proteinFrameDB::writeDBFile(string dbPath)
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

const set<int>& proteinFrameDB::getContacts(int target_i, int res_i, bool BBint) {
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
    if (aa == SeqTools::unknownIdx())
        MstUtils::error("Must set amino acid before finding interacting residues", "alignInteractingFrames::findInteractingRes");
    allInteractingFrames.clear();
    interactionData[aa] = map<int,int>();
    map<int,int>& interactionDataForAA = interactionData[aa];
    for (int i = 0; i <= 30; i++) interactionDataForAA[i] = 0;
    for (int target_idx = 0; target_idx < db.numTargets(); target_idx++) {
        const augmentedStructure& aS = db.getTarget(target_idx);
        for (int res_i = 0; res_i < aS.residueSize(); res_i++) {
            res_t res_i_aa = SeqTools::aaToIdx(aS.getResidue(res_i).getName());
            if (res_i_aa != aa) {
                continue;
            }
            const set<int> &interacting_res = db.getContacts(target_idx,res_i);
            interactionDataForAA[interacting_res.size()]++;
            for (int res_j : interacting_res) {
                mobileFrame* frame = new mobileFrame(&aS.getResidue(res_j), target_idx, res_i, res_j, aa);
                residueFrame *rFi = db.getResidueFrame(frame->getTarget(),frame->getResI());
                Transform rFi_to_ref = TransformFactory::switchFrames(refFrame, *rFi);
                rFi_to_ref.apply(*frame);
                allInteractingFrames.push_back(frame);
            }
        }
    }
}

Structure *alignInteractingFrames::getAlignedInteractingRes(int i)
{
    if ((i < 0) || (i >= allInteractingFrames.size())) MstUtils::error("Provided value " + MstUtils::toString(i) + " is out of range: (0," + MstUtils::toString(allInteractingFrames.size() - 1), "alignInteractingFrames::getAlignedInteractingRes");
    mobileFrame* frame = allInteractingFrames[i];

    return getAlignedInteractingRes(frame);
};

Structure *alignInteractingFrames::getAlignedInteractingRes(mobileFrame* frame) {
    // get the residues from the target structure
    Residue* rI = db.getResidue(frame->getTarget(),frame->getResI());
    Residue* rJ = db.getResidue(frame->getTarget(),frame->getResJ());

    Structure *interactingResiduePair = constructStructureFromResiduePair(rI, rJ, frame);

    return changeFrameToRef(interactingResiduePair,frame);
}

void alignInteractingFrames::writeAlignedInteractingResToPDB(string pdbPath, mstreal subsampleRate)
{
    // open file
    fstream pdb_out;
    MstUtils::openFile(pdb_out, pdbPath, fstream::out, "alignInteractingFrames::writeAlignedInteractingResToPDB");

    // transform interacting residues
    for (int i = 0; i < getNumInteracting(); i++)
    {
        if (MstUtils::randUnit() > subsampleRate) continue;
        Structure *interactingResiduePair = getAlignedInteractingRes(i);

        // write to file
        pdb_out << "HEADER    " << interactingResiduePair->getName() << endl;
        interactingResiduePair->writePDB(pdb_out);

        delete interactingResiduePair;
    }
    pdb_out.close();
}

void alignInteractingFrames::writeMobileFramesToBin(frameDB* frameBin) {
    for (mobileFrame* mF : allInteractingFrames) {
        frameBin->appendFrame(mF);
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