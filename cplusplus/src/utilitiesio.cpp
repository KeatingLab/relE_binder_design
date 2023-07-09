#include "utilitiesio.h"

string structurePathToName(string path) {
    return MstSys::splitPath(path,1);
}

vector<Structure*> loadStructuresFromPaths(vector<string> structure_paths, string selChainID) {
    vector<Structure*> ret;
    for (string path : structure_paths) {
        Structure* s = new Structure(path);
        string name = structurePathToName(s->getName());
        if (selChainID != "") {
            Chain* seed_chain = s->getChainByID(selChainID);
            if (seed_chain == NULL) MstUtils::error("No seed chain ("+selChainID+") found in structure: "+MstUtils::toString(name));
            Structure* selChainOnly = new Structure(*seed_chain);
            selChainOnly->setName(name);
            ret.push_back(selChainOnly);
            delete s;
        } else {
            s->setName(name);
            ret.push_back(s);
        }
    }
    return ret;
}

vector<shared_ptr<Structure>> loadSharedPointerStructuresFromPaths(vector<string> structure_paths, string selChainID) {
    vector<Structure*> structures = loadStructuresFromPaths(structure_paths, selChainID);
    vector<shared_ptr<Structure>> sp_structures;
    for (Structure* s : structures) {
        shared_ptr<Structure> sp_s = make_shared<Structure>(*s);
        sp_structures.push_back(sp_s);
    }
    return sp_structures;
}

/* --- --- --- --- --- seedBinaryFile --- --- --- --- --- */

bool seedBinaryFile::hasNext() {
    if (!readMode) MstUtils::error("hasNext not supported in write mode","seedBinaryFile::hasNext");
    return fs.peek() != EOF;
}

Structure *seedBinaryFile::next() {
    if (!readMode) MstUtils::error("next not supported in write mode","seedBinaryFile::next");
    return readNextFileSection(false).first;
}

void seedBinaryFile::skip() {
    if (!readMode) MstUtils::error("skip not supported in write mode","seedBinaryFile::skip");
    Structure* S = readNextFileSection(false).first;
    delete S;
}

void seedBinaryFile::reset() {
    if (!readMode) MstUtils::error("reset not supported in write mode","seedBinaryFile::reset");
    fs.clear(); // this is necessary in case ifstream doesn't clear eofbit
    fs.seekg(0, fs.beg);
}

void seedBinaryFile::scanFilePositions() {
    if (!readMode) MstUtils::error("scanFilePositions not supported in write mode","seedBinaryFile::scanFilePositions");
    MstTimer timer; timer.start();
    if (!_structureNames.empty()) return;
    reset();
    cout << "Scanning file positions..." << endl;
    _structureNames.clear();
    int count = 0;
    while (fs.peek() != EOF) {
        pair<Structure*,long> next;
        next = readNextFileSection(true,false);
        Structure *S = next.first;
        long pos = next.second;
        if (pos < 0) {
            cout << pos << endl;
        }
        _filePositions[S->getName()] = pos;
        _structureNames.push_back(S->getName());
        // if (count % 10000 == 0) {
        //     // report values
        //     cout << "Structure number: " << count << endl;
        //     cout << "Number of structure names:" << _structureNames.size() << endl;
        //     cout << "Number of file positions: " << _filePositions.size() << endl;
        // }
        delete S;
        count++;
    }
    timer.stop();
    cout << "Done scanning file, took " << timer.getDuration(MstTimer::msec) / 1000.0 << " sec" << endl;
}

bool seedBinaryFile::isStructureInFile(string name) {
    return _filePositions.count(name) != 0;
}

Structure *seedBinaryFile::getStructureNamed(string name) {
    if (!readMode) MstUtils::error("getStructureNamed not supported in write mode","seedBinaryFile::getStructureNamed");
    if (_filePositions.size() == 0) {
        scanFilePositions();
    }
    if (!isStructureInFile(name)) MstUtils::error("Structure "+name+" doesn't exist","seedBinaryFile::getStructureName()");
    fs.seekg(_filePositions[name], fs.beg);
    return readNextFileSection(false).first;
}

vector<Structure*> seedBinaryFile::getStructuresNamed(vector<string> names, bool strict) {
    vector<Structure*> seeds;
    for (string name: names) {
        if (!isStructureInFile(name) && strict) MstUtils::error("Structure "+name+" doesn't exist","seedBinaryFile::getStructureName()");
        seeds.push_back(getStructureNamed(name));
    }
    return seeds;
}

vector<shared_ptr<Structure>> seedBinaryFile::getStructureSPsNamed(vector<string> names, bool strict) {
    vector<shared_ptr<Structure>> seeds;
    for (string name: names) {
        if (!isStructureInFile(name) && strict) MstUtils::error("Structure "+name+" doesn't exist","seedBinaryFile::getStructureName()");
        Structure* s = getStructureNamed(name);
        seeds.emplace_back(make_shared<Structure>(*s));
        delete s;
    }
    return seeds;
}

void seedBinaryFile::jumpToStructureIndex(int idx) {
    if (!readMode) MstUtils::error("jumpToStructureIndex not supported in write mode","seedBinaryFile::jumpToStructureIndex");
    if (_structureNames.size() == 0) {
        scanFilePositions();
    }
    if (idx < 0) {
        fs.seekg(0, fs.beg);
    } else if (idx >= _structureNames.size()) {
        fs.seekg(0, fs.end);
    } else {
        string name = _structureNames[idx];
        fs.clear();
        fs.seekg(_filePositions[name], fs.beg);
    }
}

void seedBinaryFile::appendStructure(Structure *s) {
    if (readMode) MstUtils::error("appendStructure not supported in write mode","seedBinaryFile::appendStructure");
    if (s->residueSize() <= 0) MstUtils::error("Structure must have at least one residue","seedBinaryFile::appendStructure");

    if (!structure_added) {
        structure_added = true;
        if (!_append) MstUtils::writeBin(fs,string("seedBinFile")); //write the version at the top of the file
    } else {
        MstUtils::writeBin(fs,'E'); //finish the previous section
    }
    MstUtils::writeBin(fs,'S'); //start new structure section
    s->writeData(fs);
}

void seedBinaryFile::openFileStream() {
    if (readMode)
        MstUtils::openFile(fs, _filePath, fstream::in | fstream::binary, "seedBinaryFile::openFileStream");
    else if (_append)
        MstUtils::openFile(fs, _filePath, fstream::out | fstream::binary | fstream::app, "seedBinaryFile::openFileStream");
    else MstUtils::openFile(fs, _filePath, fstream::out | fstream::binary, "seedBinaryFile::openFileStream");
}

pair<Structure*,long> seedBinaryFile::readNextFileSection(bool save_metadata, bool verbose) {
    if ((fs.tellg() == 0)) {
        string version; MstUtils::readBin(fs, version);
        if (version != "seedBinFile") MstUtils::error("Provided file is not a seed binary file");
    }
    Structure* S = new Structure();
    long pos = fs.tellg();
    char sect; string prop; mstreal real_val; int dscrt_val;
    MstUtils::readBin(fs, sect);
    if (sect != 'S') MstUtils::error("The first section should be a Structure " + _filePath, "seedBinaryFile::readFileSection()");
    S->readData(fs);
    if (verbose) cout << S->getName() << endl;
    MstUtils::readBin(fs, sect);
    if (sect!='E') MstUtils::error("Section missing terminating 'E'","seedBinaryFile::readNextFileSection");
    return pair<Structure*,long>(S,pos);
}

vector<string> seedBinaryFile::getStructureNames() {
    if (!readMode) MstUtils::error("getStructureNamed not supported in write mode","seedBinaryFile::getStructureNames");
    if (_filePositions.size() == 0) {
        scanFilePositions();
    }
    return _structureNames;
}


/* --- --- --- --- --- seedPairDistributor --- --- --- --- --- */

void seedPairDistributor::setSeedGroupAsym(vector<shared_ptr<Structure>> seeds_A, vector<shared_ptr<Structure>> seeds_B) {
    /**
     * Provide seeds for group A and B directly, or by name to retrieve from the binary file.
     * Warning: there should not be any intersection between the names of the seeds in any of the four sets;
    */
    cout << "Setting asymmetric seed groups..." << endl;
    symMode = false;
    
    // vector<shared_ptr<Structure>> seeds_fromname_A = loadSeedsFromSeedBin(seed_names_A);
    // vector<shared_ptr<Structure>> seeds_fromname_B = loadSeedsFromSeedBin(seed_names_B);
    
    // // merge the direct and from-binary seeds
    // if (checkForOverlap(seeds_fromname_A,seeds_A)) MstUtils::error("Group A seeds loaded by name and those provided directly have overlapping names","seedPairDistributor::setSeedGroupAsym");
    // if (checkForOverlap(seeds_fromname_B,seeds_B)) MstUtils::error("Group B seeds loaded by name and those provided directly have overlapping names","seedPairDistributor::setSeedGroupAsym");
    
    // sort(seeds_fromname_A.begin(),seeds_fromname_A.end());
    // sort(seeds_A.begin(),seeds_A.end());
    // merge(seeds_fromname_A.cbegin(),seeds_fromname_A.cend(),
    //         seeds_A.cbegin(),seeds_A.cend(),
    //         back_inserter(seed_group_A));
    
    // sort(seeds_fromname_B.begin(),seeds_fromname_B.end());
    // sort(seeds_B.begin(),seeds_B.end());
    // merge(seeds_fromname_B.cbegin(),seeds_fromname_B.cend(),
    //         seeds_B.cbegin(),seeds_B.cend(),
    //         back_inserter(seed_group_B));

    seed_group_A = seeds_A;
    seed_group_B = seeds_B;

    cout << "seed group A: " << seed_group_A.size() << ", seed group B: " << seed_group_B.size() << endl;
    // Check for overlap between group A and B
    if (tolerateDuplicates) {
        // In this scenario, we just remove any duplicated seeds from group B
        auto seed_group_B_copy = seed_group_B;
        seed_group_B.clear();
        set<string> seed_group_A_names;
        for (shared_ptr<Structure> seed : seed_group_A) {
            seed_group_A_names.insert(seed->getName());
        }
        int count = 0;
        for (shared_ptr<Structure> seed : seed_group_B_copy) {
            // check if there is a seed with the same name in the other set
            if (seed_group_A_names.find(seed->getName()) == seed_group_A_names.end()) {
                seed_group_B.push_back(seed);
            } else {
                count++;
                cout << "(symMode = False) Seed with name " << seed->getName() << " already in group A, dropping from group B" << endl;
            }
        }
        cout << "Removed " << count << " seeds from group B" << endl;
    } else {
        if (checkForOverlap(seed_group_A,seed_group_B)) MstUtils::error("Group A and B seeds have overlapping names","seedPairDistributor::setSeedGroupAsym");
    }
    create1DJobArray();
    cout << "(symMode = false) After deduplicating, seed group A: " << seed_group_A.size() << ", seed group B: " << seed_group_B.size() << ", and total number of pairs: " << job_array.size() << endl;
}

pair<shared_ptr<Structure>,shared_ptr<Structure>> seedPairDistributor::next() {
    // increment the position
    int i = job_array[current_index].first, j = job_array[current_index].second;
    cout << current_index << " " << i << " " << j << endl;
    pair<shared_ptr<Structure>,shared_ptr<Structure>> result;
    if ((!symMode)&(j>=seed_group_B.size())) {
        j = j - seed_group_B.size();
        result = pair<shared_ptr<Structure>,shared_ptr<Structure>>(seed_group_B[j],seed_group_A[i]);
    } else {
        result = pair<shared_ptr<Structure>,shared_ptr<Structure>>(seed_group_A[i],seed_group_B[j]);
    }
    current_index++;
    return result;
}