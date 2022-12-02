#ifndef _BRIDGESEEDS_H
#define _BRIDGESEEDS_H

#include "mstfuser.h"
#include "msttypes.h"

#include "alignframes.h"
#include "hashframes.h"
#include "residueframe.h"

class seedBridge {
/**
A bridge consists of three regions
1) the region overlapping the C-terminus of seed A
2) the unique region connecting the two seeds
3) the region overlapping the N-terminus of seed B
By convention, the N and C-terminal residues of the structure come from region 1 and 3
This makes sense, since when we are trying to connect two seeds we don't have any of the connecting residues yet
*/

    public:
        seedBridge() {};
        seedBridge(int _targetID, Residue* ntermR, Residue* ctermR);

        void setRMSD(mstreal _rmsd) {rmsd = _rmsd;}
        void setUniqueID(int _uniqueID) {uniqueID = _uniqueID;}
        int getUniqueID();

        void writeToBin(ostream& ofs);
        void readFromBin(istream& ifs);

        int getTargetID() {return targetID;}
        int getNtermResIdx() {return ntermResIdx;}
        int getCtermResIdx() {return ntermResIdx+residueLength+1;}
        CartesianPoint getTerminalResidueDistances();
        int getBridgeLength() {return residueLength;}
        mstreal getRMSD() {return rmsd;}

        static CartesianPoint findCADistances(Residue* ntermR, Residue* ctermR);
        static vector<Atom*> getBridgeTerminusFromStructure(Residue* ntermR, Residue* ctermR, int terminusLength = 3, bool terminusOnly = true);

        bool operator<(seedBridge& other) {
            return this->rmsd < other.getRMSD();
        }

    private:
        int uniqueID = -1;
        int targetID = -1;
        int ntermResIdx = -1; // N-terminal residue in the original structure
        int residueLength = -1; // the total number of residues in the connecting region
        mstreal R1CaDistance = 0.0;
        mstreal R2CaDistance = 0.0;
        mstreal R3CaDistance = 0.0;
        mstreal rmsd = -1.0;
};

class seedBridgeDB {
    public:
        seedBridgeDB() {};
        seedBridgeDB(string pathToDBFile) {readDBfromFile(pathToDBFile);}
        ~seedBridgeDB() {
            for (seedBridge* sB : allBridges) delete sB;
            if (structureDB) delete structureDB;
        }

        void loadProteinStructures(string pathToStructureDB) {
            if (structureDB) delete structureDB;
            structureDB = new augmentedStructureDB();
            structureDB->readDBFile(pathToStructureDB);
        }
        void unloadProteinStructures() {delete structureDB; structureDB = nullptr;}
        
        void buildDBfromStructures(int _maxBridgeLen = 15);
        void writeDBtoFile(string pathToDBFile);
        void readDBfromFile(string pathToDBFile);
        
        seedBridge* getBridge(int uniqueID);
        const vector<seedBridge*>& getAllBridges() {return allBridges;}
        int getMaxBridgeLength() {return maxBridgeLength;}
        int getTerminusLength() {return terminusLength;}

        vector<Atom*> getBridgeTerminusFromDB(seedBridge* sB);
        Structure getBridgeAndTerminusFromDB(seedBridge* sB);

    private:
        augmentedStructureDB* structureDB = nullptr;
        vector<seedBridge*> allBridges;

        int maxBridgeLength = -1;
        int terminusLength = 3; // number of residues on either side of the bridge
};

class findSeedBridge {
    public:
        findSeedBridge(string pathToBridgeDataDB) : bridgeData(pathToBridgeDataDB) {
            loadSeedBridgeDataIntoAPV();
            info_out = new fstream;
            MstUtils::openFile(*info_out,"RMSDMatches.csv",fstream::out);
            *info_out << "name";
            for (int i = 0; i <= bridgeData.getMaxBridgeLength(); i++) *info_out << "," << i;
            *info_out << endl;
        }
        void reportAPVBoundaries();

        ~findSeedBridge() {
            info_out->close();
            delete info_out;
        }

        void setSearchQuery(Structure* nTermSeed, Structure* cTermSeed, int S1CterminusResOffset = 0, int S2NterminusResOffset = 0);

        // The initial search is performed by comparing the distance between the CA of three terminal residues on each side
        int searchSeedsByCADistance(mstreal distanceCutoff);
        vector<mstreal> getBridgeLengthDist();
        vector<seedBridge*> getBridgeInfo();

        // The matches are then verified by performing optimal superposition and calculating RMSD (this ensures proper handedness)
        void loadStructureDB(string structureDBPath) {bridgeData.loadProteinStructures(structureDBPath);}
        int verifyMatchesBySuperposition(mstreal RMSDCutoff = 0.5);
        
        // For accessing the matches
        vector<int> getVerifiedBridgeLengthDist(string name);
        vector<vector<Structure*>> getVerifiedBridgeStructures();
        // vector<Structure> getRepresentativeForEachLength();
        Structure* getAlignedBridgeStructure(seedBridge* sB);

        // For storing the matches/fused bridges
        void writeToMultiPDB(string pathToPDB, string bridgeName, int topN = -1);

        seedBridgeDB bridgeData;

    protected:
        void loadSeedBridgeDataIntoAPV();

        string getSeedBridgeName(seedBridge* sB);

    private:

        ProximitySearch PS;

        // query info
        Structure* nTermSeed = nullptr;
        Structure* cTermSeed = nullptr;
        int nTermSeedOffset = 0;
        int cTermSeedOffset = 0;
        CartesianPoint queryDistances;
        vector<Atom*> terminalResidueBBAtoms;

        vector<seedBridge*> matches;
        vector<seedBridge*> verifiedMatches;

        vector<vector<seedBridge*>> matchesByLength;

        RMSDCalculator calc;

        fstream* info_out = nullptr;
};

class fuseSeedsAndBridge {
    public:
        fuseSeedsAndBridge(int _overlapLength) : cl(false), overlapLength(_overlapLength) {
            bridge_out = new fstream;
            MstUtils::openFile(*bridge_out,"bridges.pdb",fstream::out);
            fused_out = new fstream;
            MstUtils::openFile(*fused_out,"fused.pdb",fstream::out);
        };

        ~fuseSeedsAndBridge() {
            if (!bridges.empty()) for (vector<Structure*> lengthNBridges : bridges) for (Structure* s : lengthNBridges) delete s;
            bridge_out->close();
            delete bridge_out;
            fused_out->close();
            delete fused_out;
        }

        void setClashCheck(Structure& S) {
            toCheck = S;
            clashCheck.setStructure(toCheck);
        }

        void setSeeds(Structure* _seedA, Structure* _seedB, int _seedAOffset = 0, int _seedBOffset = 0) {
            seedA =_seedA;
            seedB = _seedB;
            seedAOffset = _seedAOffset;
            seedBOffset = _seedBOffset;
        }

        void setBridgeStructures(vector<vector<Structure*>> _bridges) {
            if (!bridges.empty()) for (vector<Structure*> lengthNBridges : bridges) for (Structure* s : lengthNBridges) delete s;
            bridges = _bridges;
        }

        void writeFusedStructuresToPDB();

        vector<Structure*> clusterBridgeStructures(vector<Structure*> lengthNBridges, int Nstructures, mstreal clusterRMSD = 1.0);
        

    protected:
        vector<Residue*> getFragmentRes(Structure& seed, int startResIdx, int nRes) {
            vector<Residue*> seed_res = seed.getResidues();
            if ((startResIdx < 0)||(startResIdx>=seed_res.size())) MstUtils::error("start_residx must be in the bounds of the structure");
            if (nRes < 0) MstUtils::error("nRes must be at least 1");
            if (startResIdx+nRes>seed_res.size()) MstUtils::error("startResIdx + nRes must not be greater than the number of residues in the seed");
            return vector<Residue*>(seed_res.begin()+startResIdx,seed_res.begin()+startResIdx+nRes);
        }

        vector<int> getFragResIdx(int nTermRes, int nRes) {
            vector<int> result(nRes,0);
            iota(result.begin(),result.end(),nTermRes);
            return result;
        }

    private:
        Clusterer cl;
        Structure toCheck;
        clashChecker clashCheck;

        int overlapLength;

        Structure* seedA;
        Structure* seedB;
        int seedAOffset;
        int seedBOffset;

        vector<vector<Structure*>> bridges;

        fstream* bridge_out = nullptr;
        fstream* fused_out = nullptr;

        fusionParams params;
        fusionOutput fuserOut;

};

class potentialConnectionsSeedGraph {

};

class seedBridgeGraph {
    public:

    protected:

    private:
        int minSeedLength = 4;

};

class seedResidueBridgeGeometry {
    public:

        /**
         * @brief Computes the distance between the Ca of the residues (checks if residue is null)
         * 
         * @param Ri The first residue
         * @param Rj The second residue
         * @return mstreal 
         */
        static mstreal residueDistance(Residue* Ri, Residue* Rj) {
            if (!Ri || !Rj) MstUtils::error("Ri or Rj is null","seedResidueBridgeGeometry::residueDistance");
            Atom* RiCa = Ri->findAtom("CA");
            return  Ri->findAtom("CA")->distance(Rj->findAtom("CA"));
        };

        /**
         * @brief Finds the angle between A) the vector from the Ca of the ctermR to the ntermR and B) the "x-axis" of the local coordinate frame defined around the C-terminal residue
         * 
         * @param ctermR 
         * @param ntermR 
         * @return mstreal The cosine of the angle between A and B vectors
         */
        static mstreal cTerminalResidueDeviation(Residue* ctermR, Residue* ntermR) {
            if (!ctermR || !ntermR) MstUtils::error("ctermR or ntermR is null","seedResidueBridgeGeometry::residueDistance");
            residueFrame ctermRF(ctermR);

            // Compute the normalized vector from ctermR to ntermR
            CartesianPoint ctermRCaTontermRCa =  ntermR->findAtom("CA")->getCoor() - ctermR->findAtom("CA")->getCoor();
            
            return ctermRF.getX().cosineAngle(ctermRCaTontermRCa);
        };

        /**
         * @brief Finds the angle between A) the vector from the Ca of the ctermR to the ntermR and B) the negative "x-axis" of the local coordinate frame defined around the N-terminal residue
         * 
         * @param ctermR 
         * @param ntermR 
         * @return mstreal The cosine of the angle between A and B vectors
         */
        static mstreal nTerminalResidueDeviation(Residue* ctermR, Residue* ntermR) {
            if (!ctermR || !ntermR) MstUtils::error("ctermR or ntermR is null","seedResidueBridgeGeometry::residueDistance");
            residueFrame ntermRF(ntermR);

            // Compute the normalized vector from ctermR to ntermR
            CartesianPoint ntermRCaToctermRCa =  ctermR->findAtom("CA")->getCoor() - ntermR->findAtom("CA")->getCoor();
            
            return (-ntermRF.getX()).cosineAngle(ntermRCaToctermRCa);
        };

        static CartesianPoint averageOrientationVector(Structure* S, int windowLength = 2) {
            vector<Residue*> residues = S->getResidues();
            vector<CartesianPoint> orientationVectors;
            for (int i = 0; i < residues.size() - windowLength + 1; i++) {
                Residue* R1 = residues[i]; // get first residue
                Residue* R2 = residues[i+windowLength]; // get second residue

                // find the vector between their alpha-carbons
                orientationVectors.emplace_back(R2->findAtom("CA")->getCoor() - R1->findAtom("CA")->getCoor());
            }
            // find the average of the orientation vectors
            int n = orientationVectors.size();
            CartesianPoint averageOrientation(0,0,0);
            for (CartesianPoint p : orientationVectors) {
                averageOrientation += p.getUnit();
            }
            averageOrientation /= n;
            return averageOrientation;
        }

};

#endif