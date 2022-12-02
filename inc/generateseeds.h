#ifndef _GENERATESEEDS_H
#define _GENERATESEEDS_H

#include "mstfasst.h"
#include "mstmagic.h"
#include "mstrotlib.h"
#include "mstsystem.h"
#include "msttypes.h"

#include "residuecontact.h"
#include "freesasaext.h"

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
    Structure * getStructureNamed(string name);
    vector<Structure*> getStructuresNamed(vector<string> names);
    
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

class multiEntryPDB {
    public:
        multiEntryPDB(string filePath) : _filePath(filePath) {
            MstUtils::openFile(fs,_filePath,fstream::out);
        }

        ~multiEntryPDB() {
            fs.close();
        }

        void writeToFile(const Structure &S, string name = "") {
            fs << "HEADER    ";
            if (name == "") fs << S.getName() << endl;
            else fs << name << endl;
            S.writePDB(fs);
        }

    private:
        string _filePath;
        
        fstream fs;
};

struct seedGenParams {
    int targetFlankRes = 1;
    int maxNumMatches = 1000;
    mstreal RMSDCutoff = 0.5;
    bool seqConst = true;
    int seedFlankRes = 2;
    bool writeToPDB = false;
    string contactData = "";
    mstreal minContsPerRes = 4.0; // potential contacts per seed res
};

class seedGenerator {
    public:
        seedGenerator(const Structure& target, string fasstDBPath, seedGenParams params = seedGenParams());
 
        ~seedGenerator() {
            for (auto frag : allFragmentsData) delete frag.fragment;
        }
        
        void defineBindingSiteAsSurfaceRes(mstreal rSASAthresh = 0.1);
        void setBindingSite(vector<Residue*> bindingSiteRes) {
            bindingSite = bindingSiteRes;
            cout << "Set the binding site to " << bindingSite.size() << " residues" << endl;
        }

        string generateSeeds();

        string getTargetName() {return targetName;}

        void writeBindingSiteFragments(string path);

        struct fragmentData {
            public:
                seedGenerator* parent = nullptr;
                int cenResIdxInFrag = -1;
                Residue* fragCenResInParent = nullptr;
                Structure* fragment = nullptr;
                fasstSolutionSet matches;

                string getName() {
                    if (name == "") {
                        // parentStructure_cenResChainIDNum
                        name = parent->getTargetName() + "_" + fragCenResInParent->getChainID() 
                        + MstUtils::toString(fragCenResInParent->getNum());
                    }
                    return name;
                }

            private:
                string name = "";
        };

    protected:

        void getBindingSiteFragments();

        void findStructuralMatches();

        void generateSeedsFromMatches(seedBinaryFile& seedBin, fstream& fragOut, fstream& seedOut, multiEntryPDB* pdbOut = nullptr);

        bool seedTargetClash(Structure* seed);

    private:
        const Structure& target;
        string targetName = "";
        vector<Atom*> targetBackboneAtoms = {}; // All atoms of core residues and all backbone atoms of surface residues
        ProximitySearch targetAtomsPS;
        checkVDWRadii checker;

        potentialContacts potConts;

        vector<Residue*> bindingSite = {};

        FASST F;

        seedGenParams params;
        vector<fragmentData> allFragmentsData;

        MstTimer timer;
};

#endif