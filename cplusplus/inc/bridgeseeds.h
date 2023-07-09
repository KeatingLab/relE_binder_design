#ifndef _BRIDGESEEDS_H
#define _BRIDGESEEDS_H

#include <nlohmann/json.hpp>
using json = nlohmann::json;

#include "mstfuser.h"
#include "msttypes.h"
#include "mstsequence.h"

#include "alignframes.h"
#include "hashframes.h"
#include "residueframe.h"

#include "utilitiesio.h"

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
        seedBridgeDB(string pathToDBFile, int _maxBridgeLength, mstreal _debug_fraction = -1.0) {
            debug_fraction = _debug_fraction;
            maxBridgeLength = _maxBridgeLength;
            readDBfromFile(pathToDBFile, debug_fraction);
        }
        ~seedBridgeDB() {
            allBridges.clear();
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
        void readDBfromFile(string pathToDBFile, mstreal debug_fraction);
        
        shared_ptr<seedBridge> getBridge(int uniqueID);
        const vector<shared_ptr<seedBridge>>& getAllBridges() {return allBridges;}
        int getMaxBridgeLength() {return maxBridgeLength;}
        int getTerminusLength() {return terminusLength;}

        vector<Atom*> getBridgeTerminusFromDB(seedBridge* sB);
        Structure getBridgeAndTerminusFromDB(seedBridge* sB);

    private:
        augmentedStructureDB* structureDB = nullptr;
        vector<shared_ptr<seedBridge>> allBridges;

        int maxBridgeLength = -1;
        int terminusLength = 3; // number of residues on either side of the bridge
};

class findSeedBridge {
    public:
        findSeedBridge(string pathToBridgeDataDB, string prefix, int maxBridgeLength, mstreal debug_fraction = - 1.0) : bridgeData(pathToBridgeDataDB, maxBridgeLength) {
            loadSeedBridgeDataIntoAPV();
            info_out = new fstream;
            MstUtils::openFile(*info_out,prefix+"_RMSDMatches.csv",fstream::out);
            *info_out << "name";
            for (int i = 0; i <= bridgeData.getMaxBridgeLength(); i++) *info_out << "," << i;
            *info_out << endl;
        }
        void reportAPVBoundaries();

        ~findSeedBridge() {
            info_out->close();
            delete info_out;
        }

        void setMaxSeedDistance(mstreal _maxDistance) {maxSeedDistance = _maxDistance;}
        void setSearchQuery(const shared_ptr<Structure>& nTermSeed, const shared_ptr<Structure>& cTermSeed, int S1CterminusResOffset = 0, int S2NterminusResOffset = 0);
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
        vector<vector<shared_ptr<Structure>>> getVerifiedBridgeStructures();
        vector<vector<shared_ptr<Structure>>> getVerifiedClusteredBridgeStructures(bool verbose = true);
        // vector<Structure> getRepresentativeForEachLength();
        shared_ptr<Structure> getAlignedBridgeStructure(seedBridge* sB);

        // For storing the matches/fused bridges
        void writeToMultiPDB(string pathToPDB, string bridgeName, int topN = -1);

        seedBridgeDB bridgeData;

    protected:
        void loadSeedBridgeDataIntoAPV();

        string getSeedBridgeName(seedBridge* sB);

        vector<shared_ptr<Structure>> clusterBridgeStructures(const vector<shared_ptr<Structure>>& lengthNBridges, int Nstructures, mstreal clusterRMSD);

    private:
        mstreal maxSeedDistance = 15.0;
        mstreal debug_fraction;

        ProximitySearch PS;

        // query info
        Structure* nTermSeed = nullptr;
        Structure* cTermSeed = nullptr;
        int nTermSeedOffset = 0;
        int cTermSeedOffset = 0;
        CartesianPoint queryDistances;
        vector<Atom*> terminalResidueBBAtoms;

        vector<shared_ptr<seedBridge>> matches;
        vector<shared_ptr<seedBridge>> verifiedMatches;

        vector<vector<shared_ptr<seedBridge>>> matchesByLength;

        RMSDCalculator calc;

        fstream* info_out = nullptr;
};

struct topologyElementInfo {
    topologyElementInfo(int _topologyPosition, int _nterm_trim, int _cterm_trim, string _aaSeq, bool _isSeed) {
        topologyPosition = _topologyPosition;
        nterm_trim = _nterm_trim;
        cterm_trim = _cterm_trim;
        aaSeq = _aaSeq; 
        isSeed = _isSeed;
    }
    int topologyPosition;
    int nterm_trim;
    int cterm_trim;
    string aaSeq;
    bool isSeed; // if not seed, then it is a connector fragment
};

class topologyElement {
    public:
        topologyElement(shared_ptr<Structure> _fragment, int _topologyPosition, int _nterm_trim = 0, int _cterm_trim = 0, bool _isSeed = true) {
            // topologyPosition - the position of the N-terminal residue of the fragment within the larger topology
            // nterm_trim - residues to ignore at the N-terminus
            // cterm_trim - residues to ignore at the C-terminus
            if (_fragment->chainSize() > 1) MstUtils::error("Fragments in topology can only have a single chain","topologyElement:topologyElement");
            fragment = _fragment;
            topologyPosition = _topologyPosition;
            nterm_trim = _nterm_trim;
            cterm_trim = _cterm_trim;
            aaSeq = Sequence(*_fragment).toString(); 
            isSeed = _isSeed;
            defineName();
        }

        topologyElement(multiPDBFile& fragmentDB, string fragmentName, int _topologyPosition, int _nterm_trim, int _cterm_trim, bool _isSeed) {
            Structure* fragment_p = fragmentDB.getStructure(fragmentName);
            fragment = make_shared<Structure>(*fragment_p);
            delete fragment_p;
            topologyPosition = _topologyPosition;
            nterm_trim = _nterm_trim;
            cterm_trim = _cterm_trim;
            aaSeq = Sequence(*fragment).toString(); 
            isSeed = _isSeed;
            defineName();
        }

        bool operator<(const topologyElement& other) {
            return this->topologyPosition < other.topologyPosition;
        }

        vector<Residue*> getResidues() const {
            vector<Residue*> selectedRes;
            int nTotalRes = fragment->residueSize();
            for (int i = nterm_trim; i < nTotalRes - cterm_trim; i++) {
                Residue* R = &fragment->getResidue(i);
                selectedRes.push_back(R);
            }
            return selectedRes;
        }

        vector<int> getTopologyPositions() const {
            vector<int> selectedResPos;
            int nTotalRes = fragment->residueSize();
            for (int i = 0; i < nTotalRes - nterm_trim - cterm_trim; i++) {
                // Note the logic is slightly different than getResidues().
                selectedResPos.push_back(topologyPosition + i);
            }
            return selectedResPos;
        }

        int getResidueSize() const {
            return fragment->residueSize() - nterm_trim - cterm_trim;
        }

        void reportElement() {
            cout << "name: " << name << '\n';
            cout << "topologyPosition: " << topologyPosition << '\n';
            cout << "nterm_trim: " << nterm_trim << '\n';
            cout << "cterm_trim: " << cterm_trim << '\n';
            cout << "residueLength: " << getResidueSize() << '\n';
            cout << "isSeed: " << isSeed << endl;
        }

    shared_ptr<Structure> fragment;
    string name;
    int topologyPosition;
    int nterm_trim;
    int cterm_trim;
    string aaSeq;
    bool isSeed; // if not seed, then it is a connector fragment

    private:

        void defineName() {
            // Add appropriate prefix ('seed_' or 'connector_') if it has not already been added
            string seed_prefix, conn_prefix, fragmentName;
            fragmentName = fragment->getName();
            seed_prefix = fragmentName.substr(0,5);
            conn_prefix = fragmentName.substr(0,10);
            if (seed_prefix == "seed-") {
                if (!isSeed) MstUtils::error("Fragment ("+fragmentName+") has prefix "+seed_prefix+", but has type 'connector'","topologyElement::defineName");
                name = fragmentName;
            } else if (conn_prefix == "connector-") {
                if (isSeed) MstUtils::error("Fragment ("+fragmentName+") has prefix "+conn_prefix+", but has type 'seed'","topologyElement::defineName");
                name = fragmentName;
            } else {
                // No prefix: add it
                if (isSeed) name = "seed-" + fragmentName;
                else name = "connector-" + fragmentName;
            }
        }
};

class fragmentTopology {
    public:
        fragmentTopology() {};

        fragmentTopology(vector<shared_ptr<topologyElement>> _fragments) {
            fragments = _fragments;
            updateCurrentLength();
        }

        void addMultipleTopologyElements(vector<shared_ptr<topologyElement>> new_fragments, int overlap = 3, bool cterm = true) {
            // NOTE: this method ensures that the topology remains valid after adding the new elements
            // It's different from methods that add fragments, as the topology elements already exist and we're just trying to add them
            // This method respects the internal distance in topology positions within a set of fragments (if they are shifted, they're all shifted together)
            if (new_fragments.empty()) MstUtils::error("Passed an empty vector","addMultipleTopologyElements");
            // Make sure the new fragments are sorted N to C
            sortTopology(new_fragments);
            int nterm_pos = new_fragments.front()->topologyPosition;
            if (nterm_pos != 0) MstUtils::error("If adding topologyElements, the new elements must begin at 0","fragmentTopology::addMultipleTopologyElements");

            if (currentLength == 0) {
                // Doesn't matter: add to C-terminus starting at topPos 0
                for (shared_ptr<topologyElement> tE : new_fragments) {
                    shared_ptr<topologyElement> newtE = make_shared<topologyElement>(*tE);
                    if (currentLength > 0) checkIfFragmentIsValidExtension(newtE->topologyPosition,newtE->getResidueSize());
                    fragments.push_back(newtE);
                    updateCurrentLength();
                }
            } else if (cterm) {
                // C-term extension
                int topologyPosAdjust = currentLength - overlap;
                for (shared_ptr<topologyElement> tE : new_fragments) {
                    shared_ptr<topologyElement> newtE = make_shared<topologyElement>(*tE);
                    newtE->topologyPosition += topologyPosAdjust;
                    checkIfFragmentIsValidExtension(newtE->topologyPosition,newtE->getResidueSize());
                    fragments.push_back(newtE);
                    updateCurrentLength();
                }
            } else {
                // N-term extension
                vector<shared_ptr<topologyElement>> new_fragments_copy;
                for (shared_ptr<topologyElement> tE : new_fragments) new_fragments_copy.push_back(make_shared<topologyElement>(*tE));

                // Add C-terminal fragments first...
                for (int i = 0; i < new_fragments_copy.size(); i++) {
                    // set topologyPos to the difference between the current position and next position 
                    // (e.g. remember the distance between the fragment N-terminal positions)
                    if (i == new_fragments.size()-1) new_fragments[i]->topologyPosition = overlap - new_fragments[i]->getResidueSize();
                    new_fragments[i]->topologyPosition = new_fragments[i]->topologyPosition - new_fragments[i+1]->topologyPosition;
                }
                int topologyPosAdjust = currentLength - overlap;
                for (auto it = new_fragments_copy.rbegin(); it != new_fragments_copy.rend(); it++) {
                    shared_ptr<topologyElement> newtE = *it;
                    checkIfFragmentIsValidExtension(newtE->topologyPosition,newtE->getResidueSize());
                    for (shared_ptr<topologyElement> tE : fragments) tE->topologyPosition = tE->topologyPosition - newtE->topologyPosition;
                    fragments.push_back(newtE);
                    updateCurrentLength();
                }
            }
        }

        void addInitialFragment(shared_ptr<Structure> s, int overlap = 3, int nterm_trim = 0, int cterm_trim = 0) {
            if (s->residueSize() - (nterm_trim + cterm_trim) < overlap) MstUtils::error("Fragment is too short","fragmentTopology::addInitialFragment");
            addFragmentAtPos(s,0,nterm_trim,cterm_trim,true);
        }

        void addFragmentToEnd(shared_ptr<Structure> s, int overlap = 3, bool cterm = true, int nterm_trim = 0, int cterm_trim = 0, bool isSeed = true) {
            if (currentLength == 0) addFragmentAtPos(s,0,nterm_trim,cterm_trim,true);
            int topologyPos;
            if (cterm) topologyPos = currentLength - overlap;
            else topologyPos = overlap - s->residueSize(); //will induce a shift
            addFragmentAtPos(s,topologyPos,nterm_trim,cterm_trim,isSeed);
        }

        // topologyPos is relative to the current topology
        void addFragmentAtPos(shared_ptr<Structure> s, int topologyPos, int nterm_trim = 0, int cterm_trim = 0, bool isSeed = true) {
            if (currentLength == 0) {
                shared_ptr<topologyElement> tE = make_shared<topologyElement>(s,topologyPos,nterm_trim,cterm_trim,isSeed);
                fragments.push_back(tE);
            } else {
                checkIfFragmentIsValidExtension(topologyPos,s->residueSize()-cterm_trim-nterm_trim);
                if (topologyPos < 0) {
                    // N-term addition: need to shift the positions of fragments already in the topology
                    int topology_shift = -topologyPos;
                    for (shared_ptr<topologyElement> tE : fragments) tE->topologyPosition = tE->topologyPosition + topology_shift;
                    shared_ptr<topologyElement> tE = make_shared<topologyElement>(s,0,nterm_trim,cterm_trim,isSeed);
                    fragments.push_back(tE); //will sort later...
                } else {
                    // C-term addition
                    shared_ptr<topologyElement> tE = make_shared<topologyElement>(s,topologyPos,nterm_trim,cterm_trim,isSeed);
                    fragments.push_back(tE);
                }
            }
            updateCurrentLength();
        }

        bool checkIfFragmentIsInTopo(string fragmentName) {
            for (shared_ptr<topologyElement> tE : fragments) if (tE->fragment->getName() == fragmentName) return true;
            return false;
        }

        string getName() {
            sortTopology(fragments);
            string name = "";
            for (shared_ptr<topologyElement> tE : fragments) {
                name = name + tE->name;
                if (tE != fragments.back()) name = name + "_";
            }
            return name;
        }

        void reportTopology() {
            sortTopology(fragments);
            cout << "topology with name: " << getName() << endl;
            for (shared_ptr<topologyElement> tE : fragments) {
                tE->reportElement();
            }
        }

        shared_ptr<Structure> fuseTopology() {
            fusionTopology topology(currentLength);
            for (shared_ptr<topologyElement> tE: fragments) {
                vector<Residue*> residues = tE->getResidues();
                topology.addFragment(residues,tE->getTopologyPositions());
            }
            shared_ptr<Structure> fusedS = make_shared<Structure>(Fuser::fuse(topology,fOut,fParams));
            fusedS->setName(getName());
            fusedS->getChain(0).setID("0");
            return fusedS;
        }

        vector<shared_ptr<topologyElement>> getFragments() const {
            return fragments;
        }

        bool checkIfFragmentIsValidExtension(int topologyPos, int resLength) {
            // To be valid, must have at least minOverlap, but also must extend topology by at least one residue
            if (topologyPos < 0) {
                // N-term extension
                if (topologyPos + resLength < minOverlap) MstUtils::error("topology position "+MstUtils::toString(topologyPos)+" for fragment with "+MstUtils::toString(resLength)+" residues does not have sufficient overlap","fragmentTopology::checkIfFragmentIsValidExtension");
            } else {
                // C-term extension
                // has minOverlap?
                if (currentLength - topologyPos < minOverlap) MstUtils::error("topology position "+MstUtils::toString(topologyPos)+" in topology with "+MstUtils::toString(currentLength)+" residues has insufficient overlap","fragmentTopology::checkIfFragmentIsValidExtension");
                // extended by at least one residue?
                if (resLength - (currentLength - topologyPos) <= 0) MstUtils::error("topology position "+MstUtils::toString(topologyPos)+" for fragment with "+MstUtils::toString(resLength)+" residues in topology with "+MstUtils::toString(currentLength)+" residues does not increase the length of the topology","fragmentTopology::checkIfFragmentIsValidExtension");
            }
            return true;
        }

        void addFragmentsInTopoToDB(multiPDBFile& file) {
            for (shared_ptr<topologyElement>& tE: fragments) {
                file.addStructure(*tE->fragment,tE->name);
            }
        }

    protected:
        void updateCurrentLength() {
            sortTopology(fragments);
            currentLength = fragments.back()->topologyPosition + fragments.back()->getResidueSize();
        }

        void sortTopology(vector<shared_ptr<topologyElement>>& _fragments) {
            std:sort(_fragments.begin(),_fragments.end(),
                [](const shared_ptr<topologyElement>& lhs, const shared_ptr<topologyElement>& rhs) 
                {
                    return *(lhs.get()) < *(rhs.get());
                });
        }

    private:
        vector<shared_ptr<topologyElement>> fragments;
        int currentLength = 0;

        int minOverlap = 3;

        fusionOutput fOut;
        fusionParams fParams;
};

class topologyDB {
    public:
        topologyDB() {;}

        void addTopologyToDB(fragmentTopology& fT) {
            // Add all fragments corresponding to DB
            for (shared_ptr<topologyElement> tE : fT.getFragments()) {
                data[fT.getName()][tE->name] = {{"topologyPosition",tE->topologyPosition},
                                                {"nterm_trim",tE->nterm_trim},
                                                {"cterm_trim",tE->cterm_trim},
                                                {"aaSeq",tE->aaSeq},
                                                {"isSeed",tE->isSeed}};
            }
        }

        void readDBFromJSON(string path, bool append = true) {
            if (append) {
                fstream json_in;
                json new_data;
                MstUtils::openFile(json_in,path,fstream::in);
                json_in >> new_data;
                data.update(new_data);
            } else {
                fstream json_in;
                MstUtils::openFile(json_in,path,fstream::in);
                json_in >> data;
            }
        }

        void writeDBtoJSON(string path) {
            cout << "Writing JSON to " << path << "..." << endl;
            fstream json_out;
            MstUtils::openFile(json_out,path,fstream::out);
            json_out << data;
        }

        bool isTopoInDB(string topologyName) {
            if (data.find(topologyName) != data.end()) return true;
            return false;
        }

        shared_ptr<fragmentTopology> getTopologyFromDB(string topologyName, multiPDBFile& fragmentDB) {
            if (data.find(topologyName) == data.end()) MstUtils::error("topology not in DB","topologyDB::getTopologyFromDB");
            shared_ptr<fragmentTopology> fT = make_shared<fragmentTopology>();
            vector<shared_ptr<topologyElement>> fragments;
            // iterate over fragments
            for (json::iterator it = data[topologyName].begin(); it != data[topologyName].end(); it++) {
                int topologyPosition = it.value()["topologyPosition"];
                int nterm_trim = it.value()["nterm_trim"];
                int cterm_trim = it.value()["cterm_trim"];
                string aaSeq = it.value()["aaSeq"];
                bool isSeed = it.value()["isSeed"];
                shared_ptr<topologyElement> tE = make_shared<topologyElement>(fragmentDB,it.key(),topologyPosition,nterm_trim,cterm_trim,isSeed);
                // tE->reportElement();
                fragments.push_back(tE);
            }
            // cout << "Loaded " << fragments.size() << " fragments" << endl;
            fT->addMultipleTopologyElements(fragments);
            return fT;
        }

    private:
        json data;
};

class fuseSeedsAndBridge {
    public:
        fuseSeedsAndBridge(int _overlapLength, string prefix) : cl(false), overlapLength(_overlapLength) {
            bridge_out = new fstream;
            
            MstUtils::openFile(*bridge_out,prefix+"_bridges.pdb",fstream::out);
            fused_out = new fstream;
            MstUtils::openFile(*fused_out,prefix+"_fused.pdb",fstream::out);
        };

        ~fuseSeedsAndBridge() {
            if (!bridges.empty()) bridges.clear();
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

        void setBridgeStructures(const vector<vector<shared_ptr<Structure>>>& _bridges) {
            if (!bridges.empty()) bridges.clear();
            bridges = _bridges;
        }

        vector<shared_ptr<Structure>> fuse();

        void writeFusedStructuresToPDB();

        vector<shared_ptr<Structure>> clusterBridgeStructures(const vector<shared_ptr<Structure>>& lengthNBridges, int Nstructures, mstreal clusterRMSD = 1.0);

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
        int maxTerminiDistance = 15;

        Structure* seedA;
        Structure* seedB;
        int seedAOffset;
        int seedBOffset;

        vector<vector<shared_ptr<Structure>>> bridges;

        fstream* bridge_out = nullptr;
        fstream* fused_out = nullptr;

        fusionParams params;
        fusionOutput fuserOut;

        vector<shared_ptr<Structure>> all_fused;
};

// class potentialConnectionsSeedGraph {

// };

// class seedBridgeGraph {
//     public:

//     protected:

//     private:
//         int minSeedLength = 4;

// };

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

        static CartesianPoint averageOrientationVector(Structure* S, int windowLength = 4) {
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