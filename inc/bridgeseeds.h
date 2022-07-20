#ifndef _BRIDGESEEDS_H
#define _BRIDGESEEDS_H

#include "msttypes.h"

#include "alignframes.h"
#include "hashframes.h"
#include "residueframe.h"

class seedBridge {
    public:
        seedBridge() {};
        seedBridge(int _targetID, Residue* ntermR, Residue* ctermR);

        void setUniqueID(int _uniqueID) {uniqueID = _uniqueID;}
        int getUniqueID();

        void writeToBin(ostream& ofs);
        void readFromBin(istream& ifs);

        int getTargetID() {return targetID;}
        int getNtermResIdx() {return ntermResIdx;}
        int getCtermResIdx() {return ntermResIdx+residueLength+1;}
        CartesianPoint getTerminalResidueDistances();
        int getBridgeLength() {return residueLength;}

        static CartesianPoint findCADistances(Residue* ntermR, Residue* ctermR);
        static vector<Atom*> getBridgeTerminusFromStructure(Residue* ntermR, Residue* ctermR, int terminusLength = 3, bool terminusOnly = true);

    private:
        int uniqueID = -1;
        int targetID = -1;
        int ntermResIdx = -1;
        int residueLength = -1; // the total number of residues in the bridge spanning the N/C-terminal residues (which are not included)
        mstreal R1CaDistance = 0.0;
        mstreal R2CaDistance = 0.0;
        mstreal R3CaDistance = 0.0;
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
            structureDB = new proteinFrameDB();
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
        proteinFrameDB* structureDB = nullptr;
        vector<seedBridge*> allBridges;

        int maxBridgeLength = -1;
        int terminusLength = 3; // number of residues on either side of the bridge
};

class findSeedBridge {
    public:
        findSeedBridge(string pathToBridgeDataDB) : bridgeData(pathToBridgeDataDB) {
            loadSeedBridgeDataIntoAPV();
        }
        void reportAPVBoundaries();

        void setSearchQuery(Structure* S1, Structure* S2);

        // The initial search is performed by comparing the distance between the CA of three terminal residues on each side
        void searchSeedsByCADistance(mstreal distanceCutoff);
        vector<mstreal> getBridgeLengthDist();
        vector<seedBridge*> getBridgeInfo();

        // The matches are then verified by performing optimal superposition and calculating RMSD (this ensures proper handedness)
        void loadStructureDB(string structureDBPath) {bridgeData.loadProteinStructures(structureDBPath);}
        void verifyMatchesBySuperposition(mstreal RMSDCutoff = 0.5);
        vector<mstreal> getVerifiedBridgeLengthDist();
        vector<Structure> getVerifiedBridgeStructures(int bridgeLength = -1);

    protected:
        void loadSeedBridgeDataIntoAPV();

    private:
        seedBridgeDB bridgeData;

        ProximitySearch PS;

        // query info
        Structure* nTermSeed = nullptr;
        Structure* cTermSeed = nullptr;
        CartesianPoint queryDistances;
        vector<Atom*> terminalResidueBBAtoms;

        vector<seedBridge*> matches;
        vector<seedBridge*> verifiedMatches;

        RMSDCalculator calc;
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