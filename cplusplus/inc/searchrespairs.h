#ifndef _SEARCHRESPAIRS_H
#define _SEARCHRESPAIRS_H

#include "msttypes.h"

#include "alignframes.h"
#include "hashframes.h"
#include "residueframe.h"

class distance3AngleTable {
    public:
        distance3AngleTable(mstreal _maxVal, mstreal _dCut, mstreal _angleCut) : maxVal(_maxVal), dCut(_dCut), angleCut(_angleCut) {
            nBins = ceil(2*((maxVal)/dCut));
            cout << "nBins distance: " << nBins << endl;
            binLength = maxVal/mstreal(nBins-1); // should be an extra bin at the end
            cout << "binLength: " << binLength << endl;
            for (int i = 0; i < nBins; i++) {
                // cout << "angleCut: " << angleCut << endl;
                cout << "nBins angles: " << 2*ceil(180/angleCut) << endl;
                bins.emplace_back(0,0,0,180,180,180,ceil(2*180/angleCut));
            }
        }

        void addResPairToTable(resPair* rP, int tag) {
            if (rP->getCaDistance() > 0.0) {
                // int idx = val2Idx(rP->getCaDistance());
                // cout << "bin idx: " << val2Idx(rP->getCaDistance()) << endl;
                // cout << bins[idx].getXHigh() << " " << bins[idx].getYHigh() << " " << bins[idx].getZHigh() << endl;
                bins[val2Idx(rP->getCaDistance())].addPoint(rP->getAngles(),tag);
                idx2distance[tag] = rP->getCaDistance();
            }
        }

        vector<int> searchResPair(resPair* rP);

        vector<int> getBinIdxRange(mstreal val, mstreal tol) {
            int lowerBin = val2Idx(val-tol);
            int higherBin = val2Idx(val+tol);
            cout << "lowerBin: " << lowerBin << " higherBin: " << higherBin << endl;
            vector<int> range;
            for (int i = lowerBin; i <= higherBin; i++) range.push_back(i);
            return range;
        }
    
    protected:
        int val2Idx(mstreal val) {
            if (bins.empty()) MstUtils::error("Bins are empty","lookup1D::val2Idx");
            if (val > maxVal) MstUtils::error("Provided value: "+MstUtils::toString(val)+" greater than the maximum value in the table","lookup1D::val2Idx");
            return floor(val/binLength);
        }

    private:
        mstreal maxVal;
        int nBins;
        mstreal binLength;
        vector<ProximitySearch> bins;
        map<int,mstreal> idx2distance;

        mstreal dCut;
        mstreal angleCut;
};

class findResPairs {
    public:
        findResPairs(string resPairDBPath, mstreal _dCut = 0.25, mstreal _angleCut = 30, mstreal _rmsdCut = 0.25);

        ~findResPairs() {
            for (resPair* rP : allResPairs) delete rP;
            // if (resPairMap != nullptr) delete resPairMap;
        }

        void setQuery(Residue* Ri, Residue* Rj) {
            queryRP = resPair(Ri,Rj);
        }

        resPair& getQuery() {return queryRP;}

        int searchForMatches(bool verbose = true);

        resPair* getMatch(int matchIdx);
        vector<resPair*> getMatches() {return verifiedMatches;}
        int getNumMatchesWithResidueType(bool Ri = true);
        vector<int> getNumMatchesByAAType(bool Ri = true);

    protected:

    private:
        resPairDB DB;
        vector<resPair*> allResPairs;
        
        // distance3AngleTable* resPairMap = nullptr;
        ProximitySearch PS;

        resPair queryRP;
        mstreal dCut;
        mstreal angleCut;
        mstreal rmsdCut;

        vector<resPair*> matches;
        vector<resPair*> verifiedMatches;

        RMSDCalculator calc;
};

class findkNNResPairs {
    public:
        findkNNResPairs(string resPairDBPath, int _k = 30, bool _noSearch = false) : resPairSearcher(resPairDBPath,0.75,30,0.25), k_max(_k), noSearch(_noSearch) {;}

        void setStructure(Structure* _S) {
            k = k_max; // reset to original value
            if (_S->residueSize() < k) {
                if (_S->residueSize() < 2) MstUtils::error("Structure must have at least 2 residues");
                k = _S->residueSize();
                cout << "Warning: structure has less than 30 residues, k will temporary be decreased to " << k << endl;
                // MstUtils::error("Structure must have at least 30 residues");
            }
            S = _S;
            resCa = AtomPointerVector(MiscTools::getBackboneCaAtoms(S->getResidues()));
            PS = ProximitySearch(resCa,5.0);
            kNN.clear();
            kNN_to_matches.clear();
        }

        set<pair<Residue*,Residue*>> findkNN() {
            set<pair<Residue*,Residue*>> allNeighbors;
            for (Residue* Ri : S->getResidues()) {
                vector<Residue*> NN = findkNN(Ri);
                for (Residue* Rj : NN) {
                    if (Ri == Rj) continue;
                    pair<Residue*,Residue*> res_pair = Ri->getResidueIndex() <= Rj->getResidueIndex() ? pair<Residue*,Residue*>(Ri,Rj) : pair<Residue*,Residue*>(Rj,Ri);
                    allNeighbors.insert(res_pair);
                }
            }
            kNN = allNeighbors;
            cout << "Found " << kNN.size() << " neighbors" << endl;
            return allNeighbors;
        }

        void searchkNN() {
            if (kNN.empty()) MstUtils::error("Cannot search for matches until neighbors have been defined");
            for (pair<Residue*,Residue*> pair : kNN) {
                Residue* Ri = pair.first;
                Residue* Rj = pair.second;

                resPairSearcher.setQuery(Ri,Rj);
                int n_matches;
                if (noSearch) {
                    n_matches = 0;
                } else {
                    n_matches = resPairSearcher.searchForMatches();
                }
                kNN_to_matches[pair] = n_matches;
            }
        }

        void writeToFile() {
            fstream out;
            MstUtils::openFile(out,S->getName()+"_respair_matches.csv",fstream::out);
            out << "Ri_resIdx,Ri_chainID,Ri_resnum,Rj_resIdx,Rj_chainID,Rj_resnum,Ca_distance,n_matches" << endl;
            for (pair<Residue*,Residue*> pair : kNN) {
                Residue* Ri = pair.first;
                Residue* Rj = pair.second;

                out << Ri->getResidueIndex() << "," << Ri->getChainID() << "," << Ri->getNum() << ",";
                out << Rj->getResidueIndex() << "," << Rj->getChainID() << "," << Rj->getNum() << ",";
                out << Ri->findAtom("CA")->getCoor().distance(Rj->findAtom("CA")->getCoor()) << ",";
                out << kNN_to_matches[pair] << endl;
            }
        }

    protected:
        vector<Residue*> findkNN(Residue* R) {
            vector<Residue*> NN;

            // Search until we find 30 or more neighbors
            Atom* Ca = R->findAtom("CA");
            CartesianPoint Ca_coord = Ca->getCoor(); 
            vector<int> neighbors;
            mstreal radius = 10.0;
            while (neighbors.size() < k) {
                neighbors = PS.getPointsWithin(Ca_coord,0,radius);
                radius+=5.0;
            }

            // Sort by distance
            std::sort(neighbors.begin(),neighbors.end(),
                [this,Ca_coord](int i, int j)
            {
                return this->resCa[i]->getCoor().distance(Ca_coord) < this->resCa[j]->getCoor().distance(Ca_coord);
            });

            for (int atom : neighbors) {
                Residue* neighbor_residue = resCa[atom]->getParent();
                // cout << neighbor_residue->getChainID() << neighbor_residue->getNum() << endl;
                // cout << "distance: " << resCa[atom]->getCoor().distance(Ca_coord) << endl;
                NN.push_back(neighbor_residue);
                if (NN.size() >= k) break;
            }
            return NN;
        }

    private:
        int k_max = 30; //includes self
        int k = k_max;
        bool noSearch = false;

        findResPairs resPairSearcher;

        Structure* S = nullptr;
        ProximitySearch PS;
        AtomPointerVector resCa;

        set<pair<Residue*,Residue*>> kNN; //does not include duplicates. residues sorted: res i < res j
        map<pair<Residue*,Residue*>,int> kNN_to_matches;
};

#endif