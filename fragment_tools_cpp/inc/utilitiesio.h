#ifndef _UTILITIESIO_H
#define _UTILITIESIO_H

#include <memory>

#include "mstrotlib.h"
#include "mstsequence.h"
#include "mstsystem.h"

string structurePathToName(string path);

vector<Structure*> loadStructuresFromPaths(vector<string> structure_paths, string selChainID = "");

vector<shared_ptr<Structure>> loadSharedPointerStructuresFromPaths(vector<string> structure_paths, string selChainID = "");

/** 
 * Interface for loading seeds from file. Not responsible for managing memory after seeds are loaded
*/
class seedBinaryFile {
public:
    seedBinaryFile(string filePath, bool read = true, bool append = false): _filePath(filePath), readMode(read), _append(append){
        cout << "read mode: " << readMode << "\tappend: " << _append << "\topening file: " << filePath << endl;
        openFileStream();
    };

    seedBinaryFile(const seedBinaryFile &other): readMode(true), _filePath(other._filePath), _filePositions(other._filePositions), _structureNames(other._structureNames) {
        if (!other.readMode) MstUtils::error("Copying write-only binary file not supported","seedBinaryFile::seedBinaryFile");
        cout << "Opening file stream for copy, with " << _structureNames.size() << " loaded structure names" << endl;
        openFileStream();
    }

    ~seedBinaryFile() {
        if ((!readMode) && (structure_added)) MstUtils::writeBin(fs,'E');
        fs.close();
    }

    Structure *next();
    bool hasNext();
    void skip();
    void reset();
    bool isStructureInFile(string name);
    Structure * getStructureNamed(string name);
    vector<Structure*> getStructuresNamed(vector<string> names, bool strict = true);
    vector<shared_ptr<Structure>> getStructureSPsNamed(vector<string> names, bool strict = true);
    
    /*
     Finds the position of each structure section in the file
     */
    void scanFilePositions();

    size_t structureCount() {
        if (_filePositions.empty())
            scanFilePositions();
            reset();
        return _filePositions.size();
    }

    void jumpToStructureIndex(int idx);

    void insertStructureNames(vector<string> &names) {
        if (_filePositions.empty())
            scanFilePositions();
        names.insert(names.end(), _structureNames.begin(), _structureNames.end());
    }
    
    void appendStructure(Structure *s);
    
    vector<string> getStructureNames();

protected:
    void openFileStream();

    pair<Structure*,long> readNextFileSection(bool save_metadata = false, bool verbose = false);

private:
    string _filePath;
    bool readMode;
    bool _append;
    bool structure_added = false;
    
    fstream fs;
    long fileStart;
    unordered_map<string, long> _filePositions; // Positions of each structure in the file
    vector<string> _structureNames;
};

// class seedOneToManyDistributor {
//     // Note: object takes ownership of all seeds that are passed to it
//     // in contrast to seedPairDistributor
//     public:
//         seedOneToManyDistributor(int _workerID = 0, int _nWorkers = 1, string seedbinpath = "") {
//             workerID = _workerID;
//             nWorkers = _nWorkers;
//             if (seedbinpath != "") {
//                 seedbin = new seedBinaryFile(seedbinpath);
//                 seedbin->scanFilePositions();
//                 cout << seedbin->structureCount() << " structures in binary file" << endl;
//             }
//         }

//         ~seedOneToManyDistributor() {
//             if (seedbin != nullptr) delete seedbin;
//         };

//         // Directly provide the seed structures with the following three methods
//         void setSingleSeed(Structure* seed) {
//             single_seed = make_shared<Structure>(*seed);
//             delete seed;
//         }

//         void setSeedGroup(vector<string> seed_names) {
//             vector<Structure*> seeds = seedbin->getStructuresNamed(seed_names);
//             setSeedGroup(seeds);
//         }

//         void setSeedGroup(vector<Structure*> seeds) {
//             if (!seed_group.empty()) MstUtils::error("Seed group already set!","seedOneToManyDistributor::setSeedGroup");
//             int seed = 42;
//             shuffle(seeds.begin(),seeds.end(),std::default_random_engine(seed));
//             for (Structure* seed : seeds) {
//                 shared_ptr<Structure> seed_sp = make_shared<Structure>(*seed);
//                 seed_group.push_back(seed_sp);
//             }
//             // Add a second time for search connections in the other direction
//             for (Structure* seed : seeds) {
//                 shared_ptr<Structure> seed_sp = make_shared<Structure>(*seed);
//                 seed_group.push_back(seed_sp);
//                 delete seed;
//             }
//         };

//     bool hasNext() {
//         // Batch is done when we've iterated over all 
//         return ((current_index<max_index));
//     }

//     pair<shared_ptr<Structure>,shared_ptr<Structure>> next() {
//         int batch_size = batchSize();
//         if (current_index == -1) {
//             // Start the index according to the batch id
//             int batch_size = batchSize();
//             current_index = int(workerID * batch_size / seed_group.size());
//             min_index = current_index;
//             max_index = current_index + batch_size;
//         }
//         // increment the position
//         current_index++;
//         if ((current_index / 2) < seed_group.size() / 2) {
//             // single seed + group[i] seed
//             pair<shared_ptr<Structure>,shared_ptr<Structure>> result(single_seed,seed_group[current_index]);
//             return result;
//         } else {
//             // group[i] seed + single seed
//             pair<shared_ptr<Structure>,shared_ptr<Structure>> result(seed_group[current_index],single_seed);
//             return result;
//         }
//     }

//     int batchSize() {
//         return ceil(mstreal(seed_group.size()) / mstreal(nWorkers));
//     }

//     private:
//         int workerID = 0; // workers are 0-indexed
//         int nWorkers = 1;

//         int current_index = -1;

//         // range is [min,max), size of range depends on number of seeds in each group and number of workers
//         // idx is 0-indexed of course (to match the vector)
//         int min_index = 0;
//         int max_index = 0;

//         seedBinaryFile* seedbin = nullptr;
//         shared_ptr<Structure> single_seed;
//         vector<shared_ptr<Structure>> seed_group;
// };

class seedPairDistributor {
    // Note: object takes ownership of all seeds that are passed to it
    // in contrast to seedPairDistributor
    // TWO MODES: symmetric and asymmetric
    public:
        seedPairDistributor(int _workerID = 0, int _nWorkers = 1, bool _shuffle = false) {
            workerID = _workerID;
            nWorkers = _nWorkers;
            shuffle = _shuffle;
            // if (seedbinpath != "") {
            //     seedbin = new seedBinaryFile(seedbinpath);
            //     seedbin->scanFilePositions();
            //     cout << seedbin->structureCount() << " structures in binary file" << endl;
            // }
        }

        /**
         *  The logic for generating seed pairs depends on the composition of the seed groups:
         *  If seed group A and B are identical, it is the "symmetric" mode.
         *  If seed group A and B are non-overlapping, it is the "asymmetric" mode
        */
        // void setSeedGroupSym(vector<string> seed_names) {
        //     vector<shared_ptr<Structure>> seeds_sp = loadSeedsFromSeedBin(seed_names);
        //     setSeedGroupSym(seeds_sp);
        // }

        void setSeedGroupSym(vector<shared_ptr<Structure>> seeds) {
            cout << "Setting symmetric seed groups..." << endl;
            symMode = true;
            seed_group_A = seeds;
            seed_group_B = seeds;

            create1DJobArray();
            updateName2Seed();
            cout << "(symMode = true) seed group A: " << seed_group_A.size() << ", seed group B: " << seed_group_B.size() << ", and total number of pairs: " << job_array.size() << endl;
        }

        void setSeedGroupAsym(vector<shared_ptr<Structure>> seeds_A, vector<shared_ptr<Structure>> seeds_B);

        bool hasNext() {
            // Batch is done when we've iterated over all 
            return ((current_index<max_index));
        }

        pair<shared_ptr<Structure>,shared_ptr<Structure>> next();

        int batchSize() {
            return ceil(mstreal(job_array.size()) / mstreal(nWorkers));
        }

        vector<string> getAllSeedNames() {
            vector<string> ret;
            for (shared_ptr<Structure> sA : seed_group_A) ret.push_back(sA->getName());
            for (shared_ptr<Structure> sB : seed_group_B) ret.push_back(sB->getName());
            return ret;
        }

        shared_ptr<Structure> getSeedByName(string name) {
            if (name2seed.count(name) == 0) MstUtils::error("That seed aint done been loaded, partner","seedPairDistributor::getSeedByName");
            return name2seed[name];
        }

    protected:
        // vector<shared_ptr<Structure>> loadSeedsFromSeedBin(vector<string> seed_names) {
        //     // silly hack because I don't want to rework the seedBinaryFile class
        //     vector<Structure*> seeds = seedbin->getStructuresNamed(seed_names);
        //     vector<shared_ptr<Structure>> seeds_sp;
        //     for (Structure* s : seeds) {
        //         shared_ptr<Structure> seed_sp = make_shared<Structure>(*s);
        //         delete s;
        //         seeds_sp.push_back(seed_sp);
        //     }
        //     return seeds_sp;
        // }

        bool checkForOverlap(const vector<shared_ptr<Structure>>& seedsA, const vector<shared_ptr<Structure>>& seedsB) {
            vector<string> seed_A_names, seed_B_names, name_intersection;
            for (const shared_ptr<Structure>& s : seedsA) seed_A_names.push_back(s->getName());
            for (const shared_ptr<Structure>& s : seedsB) seed_B_names.push_back(s->getName());
            sort(seed_A_names.begin(),seed_A_names.end());
            sort(seed_B_names.begin(),seed_B_names.end());
            set_intersection(seed_A_names.begin(), seed_A_names.end(), seed_B_names.begin(), seed_B_names.end(),
                          std::back_inserter(name_intersection));
            if (name_intersection.size() > 0) return true;
            return false;
        }

        void create1DJobArray() {
            job_array.clear();
            if (symMode) {
                for (int i = 0; i < seed_group_A.size(); i++) {
                    for (int j = 0; j < seed_group_B.size(); j++) {
                        job_array.emplace_back(i,j);
                    }
                }
            } else {
                // Create extra pairs corresponding to the reverse order (e.g. we want A->B and B->A)
                for (int i = 0; i < seed_group_A.size(); i++) {
                    for (int j = 0; j < 2*seed_group_B.size(); j++) {
                        job_array.emplace_back(i,j);
                    }
                }
            }
            if (shuffle) {
                cout << "Shuffling the job array..." << endl; 
                auto rng = std::default_random_engine(42);
                std::shuffle(job_array.begin(),job_array.end(), rng);
            }

            int batch_size = batchSize();
            // Start the index according to the batch id
            current_index = int(workerID * batch_size);
            min_index = current_index;
            max_index = min(current_index + batch_size,int(job_array.size()));
            cout << "nWorkers: " << nWorkers << ", workerID: " << workerID << ", current_index: " << current_index << ", max_index: " << max_index << endl;
        }

        void updateName2Seed() {
            for (shared_ptr<Structure> sA : seed_group_A) name2seed[sA->getName()] = sA;
            for (shared_ptr<Structure> sB : seed_group_B) name2seed[sB->getName()] = sB;
        }

    private:
        int workerID = 0; // workers are 0-indexed
        int nWorkers = 1;
        bool shuffle = false;
        bool tolerateDuplicates = true;
        bool symMode = true;

        int current_index = -1;

        // range is [min,max), size of range depends on number of seeds in each group and number of workers
        // idx is 0-indexed of course (to match the vector)
        int min_index = 0;
        int max_index = 0;

        // seedBinaryFile* seedbin = nullptr;
        vector<shared_ptr<Structure>> seed_group_A;
        vector<shared_ptr<Structure>> seed_group_B;
        vector<pair<int,int>> job_array;
        map<string,shared_ptr<Structure>> name2seed;
};

// class seedPairDistributor {
//     // Note: object takes ownership of all seeds that are passed to it
//     public:
//         seedPairDistributor(int _workerID = 0, int _nWorkers = 1, string seedbinpath = "") {
//             workerID = _workerID;
//             nWorkers = _nWorkers;
//             if (seedbinpath != "") {
//                 seedbin = new seedBinaryFile(seedbinpath);
//                 seedbin->scanFilePositions();
//                 cout << seedbin->structureCount() << " structures in binary file" << endl;
//             }
//         }

//         ~seedPairDistributor() {
//             for (Structure* seed: all_loaded_seeds) delete seed;
//             if (seedbin != nullptr) delete seedbin;
//         };

//         // Directly provide the seed structures with the following three methods
//         void setSeeds(vector<Structure*> seeds) {
//             seed_group_a = seeds;
//             seed_group_b = seeds;
//             for (Structure* seed: seeds) all_loaded_seeds.insert(seed);
//         };

//         void setSeedsGroupA(vector<Structure*> seeds) {
//             seed_group_a = seeds;
//             for (Structure* seed: seeds) all_loaded_seeds.insert(seed);
//         }

//         void setSeedsGroupB(vector<Structure*> seeds) {
//             seed_group_b = seeds;
//             for (Structure* seed: seeds) all_loaded_seeds.insert(seed);
//         }

//         // Load the seeds from a binary file with the following methods
//         void setSeeds(vector<string> seed_names) {
//             vector<Structure*> seeds = seedbin->getStructuresNamed(seed_names);
//             setSeeds(seeds);
//         };

//         void setSeedsGroupA(vector<string> seed_names) {
//             vector<Structure*> seeds = seedbin->getStructuresNamed(seed_names);
//             setSeedsGroupA(seeds);
//         }

//         void setSeedsGroupB(vector<string> seed_names) {
//             vector<Structure*> seeds = seedbin->getStructuresNamed(seed_names);
//             setSeedsGroupB(seeds);
//         }

//     bool hasNext() {
//         // Batch is done when we've iterated over all 
//         return ((group_a_current_index<group_a_max_index)|(group_b_current_index<group_b_max_index));
//     }

//     pair<Structure*,Structure*> next() {
//         if ((group_a_current_index == -1)|(group_b_current_index == -1)) {
//             // we round up when determining the batch size, so the final batch is smaller

//             // There are A x B jobs total (imagine a grid). We use the worker ID and batch size to get our index in the total number of 
//             // jobs (which is some value less than A * B - 1). We convert this into a position in the grid, e.g., we get the row/column
//             // of A and B respectively
//             int batch_size = batchSize();
//             group_a_current_index = int(workerID * batch_size / seed_group_b.size());
//             group_b_current_index = int(workerID * batch_size % seed_group_b.size());
//             // group_a_min_index = group_a_current_index;
//             // group_b_min_index = group_b_current_index;
//             group_a_max_index = min(int(((workerID + 1) * batch_size) / seed_group_b.size()),int(seed_group_a.size()));
//             group_b_max_index = min(int(((workerID + 1) * batch_size) % seed_group_b.size()),int(seed_group_b.size()));
//         }
//         pair<Structure*,Structure*> result(seed_group_a[group_a_current_index],seed_group_b[group_b_current_index]);
//         // Now increment the position in the grid
//         group_b_current_index++;
//         if (group_b_current_index == seed_group_b.size()) {
//             group_b_current_index = 0;
//             group_a_current_index++;
//         }
//         return result;
//     }

//     int batchSize() {
//         return ceil(mstreal(seed_group_a.size() * seed_group_b.size()) / mstreal(nWorkers));
//     }

//     private:
//         int workerID = 0; // workers are 0-indexed
//         int nWorkers = 1;

//         int group_a_current_index = -1;
//         int group_b_current_index = -1;

//         // range is [min,max), size of range depends on number of seeds in each group and number of workers
//         // idx is 0-indexed of course (to match the vector)
//         int group_a_min_index = 0;
//         int group_a_max_index = 0;
//         int group_b_min_index = 0;
//         int group_b_max_index = 0;

//         seedBinaryFile* seedbin = nullptr;
//         vector<Structure*> seed_group_a;
//         vector<Structure*> seed_group_b;
//         set<Structure*> all_loaded_seeds;
        
// };

class multiPDBFile {
    public:
        multiPDBFile(string _path, bool _read_mode = true) {
            path = _path;
            fileh = new fstream;
            read_mode = _read_mode;
            if (read_mode) {
                MstUtils::openFile(*fileh,path,fstream::in);
                getFilePositions();
            } else {
                MstUtils::openFile(*fileh,path,fstream::out);
            }
        }

        ~multiPDBFile() {
            if (!read_mode) writeFileTerm();
            fileh->close();
            delete fileh;
        }

        void addStructure(const Structure& s, bool strict = false) {
            addStructure(s,s.getName(),strict);
        }

        void addStructure(const Structure& s, string structure_name, bool strict = false) {
            if (read_mode) MstUtils::error("Cannot add structure in read mode","multiPDBFile::addStructure");
            if (name2idx.find(structure_name) == name2idx.end()) {
                name2idx[structure_name] = structure_filepos.size() - 1;
                structure_filepos.push_back(addStructureToFile(s,structure_name));
            } else {
                if (strict) MstUtils::error("Tried to add: "+s.getName()+", but structure with that name already exists","multiPDBFile::addStructure");
            }
        }

        bool addStructureLines(const vector<string>& lines, bool strict = false) {
            // Note: only use this if you trust that your input is valid
            if (read_mode) MstUtils::error("Cannot add structure in read mode","multiPDBFile::addStructureLines");
            string structure_name = getNameFromHeader(lines[0]);
            if (name2idx.find(structure_name) == name2idx.end()) {
                name2idx[structure_name] = structure_filepos.size() - 1;
                structure_filepos.push_back(addStructureLinesToFile(lines));
                return true;
            } else {
                if (strict) MstUtils::error("Tried to add: "+structure_name+", but structure with that name already exists","multiPDBFile::addStructure");
                return false;
            }
        }

        int countStructures(bool verbose = false) {
            if (structure_filepos.size() == 0) getFilePositions();
            // if (verbose) cout << "Loaded multipdb with " << structure_filepos.size() << " structures" << endl;
            return structure_filepos.size();
        }

        Structure* getStructure(string name) {
            if (!read_mode) MstUtils::error("Cannot get structure in write mode","multiPDBFile::getStructure");
            if (structure_filepos.size() == 0) getFilePositions();
            if (name2idx.find(name) == name2idx.end()) MstUtils::error(name+" (name) is not in database","multPDBFile::getStructure");
            return getStructure(name2idx[name]);
        }

        Structure* getStructure(int idx) {
            if (!read_mode) MstUtils::error("Cannot get structure in write mode","multiPDBFile::getStructure");
            if (structure_filepos.size() == 0) getFilePositions();
            if ((idx < 0)||(idx > structure_filepos.size())) MstUtils::error(MstUtils::toString(idx)+" (index) is out of range","multiPDBFile::getStructure");
            // cout << idx << " " << structure_filepos[idx] << endl;
            return readStructureFromFile(structure_filepos[idx]);
        }

        shared_ptr<Structure> getStructureSP(string name) {
            Structure* s = getStructure(name);
            shared_ptr<Structure> s_sp = make_shared<Structure>(*s);
            delete s;
            return s_sp;
        }

        shared_ptr<Structure> getStructureSP(int index) {
            Structure* s = getStructure(index);
            shared_ptr<Structure> s_sp = make_shared<Structure>(*s);
            delete s;
            return s_sp;
        }

        vector<shared_ptr<Structure>> getStructuresSP(vector<string> names) {
            vector<shared_ptr<Structure>> ret;
            for (string name : names) ret.push_back(getStructureSP(name));
            return ret;
        }

        vector<shared_ptr<Structure>> loadAllStructuresSP() {
            vector<shared_ptr<Structure>> ret;
            for (int i = 0; i < countStructures(true); i++) ret.push_back(getStructureSP(i));
            return ret;
        }

        vector<string> getStructureLines(string name) {
            if (!read_mode) MstUtils::error("Cannot get structure lines in write mode","multiPDBFile::getStructureLines");
            if (structure_filepos.size() == 0) getFilePositions();
            if (name2idx.find(name) == name2idx.end()) MstUtils::error(name+" (name) is not in database","multPDBFile::getStructureLines");
            return getStructureLines(name2idx[name]);
        }

        vector<string> getStructureLines(int idx) {
            if (!read_mode) MstUtils::error("Cannot get structure lines in write mode","multiPDBFile::getStructureLines");
            if (structure_filepos.size() == 0) getFilePositions();
            if ((idx < 0)||(idx > structure_filepos.size())) MstUtils::error(MstUtils::toString(idx)+" (index) is out of range","multiPDBFile::getStructureLines");
            // cout << "idx: " << idx << " filepos: " << structure_filepos[idx] << endl;
            return readStructureLinesFromFile(structure_filepos[idx]);
        }

        void copyStructuresToNewMultiPDB(multiPDBFile& newFile) {
            if (structure_filepos.size() == 0) getFilePositions();
            int count = 0;
            for (int i = 0; i < countStructures(); i++) {
                vector<string> structure_lines = getStructureLines(i);
                if (structure_lines.empty()) break;
                bool added = newFile.addStructureLines(structure_lines);
                if (added) count++;
            }
            cout << "Added " << count << " structures to the new multiPDB" << endl;
        }

    protected:

        Structure* readStructureFromFile(unsigned long long fpos) {
            fileh->seekg(fpos);
            // cout << "fpos " << fileh->tellg() << " " << fpos << endl;
            string header;
            getline(*fileh,header);
            // cout << "header line: " << header << endl;
            string name = getNameFromHeader(header);
            Structure* s = new Structure;
            s->readPDB(*fileh);
            s->setName(name);
            return s;
        }

        vector<string> readStructureLinesFromFile(unsigned long long fpos) {
            if (fpos >= 0) fileh->seekg(fpos);
            vector<string> all_lines;
            string line;
            getline(*fileh,line);
            if ((line == "FILEEND")||(line.size() == 0)) return all_lines;
            all_lines.push_back(line);
            while (line != "END") {
                // cout << fileh->tellg() <<" line: " << line << endl;
                getline(*fileh,line);
                all_lines.push_back(line);
            }
            return all_lines;
        }

        string skipThroughPDBEntry() {
            string line;
            getline(*fileh,line);
            if ((line == "FILEEND")||(line.size() == 0)) return "";
            string name = getNameFromHeader(line);
            while (line != "END") {
                getline(*fileh,line);
            }
            return name;
        } 

        string getNameFromHeader(string header) {
            if (header.size() < 10) MstUtils::error("Line is too short to be a header: "+header,"multiPDBFile::getNameFromHeader");
            return header.substr(10);
        }

        int addStructureToFile(const Structure& s, string name) {
            int fpos = fileh->tellg();
            writeHeader(name);
            s.writePDB(*fileh);
            return fpos;
        }

        int addStructureLinesToFile(const vector<string>& lines) const {
            int fpos = fileh->tellg();
            for (const string& line : lines) *fileh << line << "\n";
            return fpos;
        }

        void writeHeader(string name) {
            *fileh << "HEADER    " << name << endl;
        }

        void writeFileTerm() {
            *fileh << "FILEEND" << endl; // line will be stripped when read by PyMol
        }

        void getFilePositions() {
            fileh->seekg(0);
            int fpos = fileh->tellg();
            string name = skipThroughPDBEntry();
            while (name != "") {
                name2idx[name] = structure_filepos.size();
                structure_filepos.push_back(fpos);
                // cout << "fpos: " << fpos << ", name: " << name << endl;
                fpos = fileh->tellg();
                if (fpos < structure_filepos.back()) MstUtils::error("File position has overflowed the unsigned long long variable");
                name = skipThroughPDBEntry();
            }
            fileh->clear();
        }

    private:
        string path;
        fstream* fileh;
        bool read_mode;

        vector<unsigned long long> structure_filepos;
        map<string,int> name2idx;
};

#endif