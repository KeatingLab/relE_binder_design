#include "searchrespairs.h"

/* --- --- --- --- --- distance3AngleTable --- --- --- --- --- */

// vector<int> distance3AngleTable::searchResPair(resPair* rP) {
//     mstreal query_distance = rP->getCaDistance();
//     vector<int> binIdxs = getBinIdxRange(query_distance,dCut);
//     vector<int> result;
//     for (int idx : binIdxs) {
//         vector<int> angle_matches = bins[idx].getPointsWithin(rP->getAngles(),0,angleCut,true);
//         cout << angle_matches.size() << " angle matches in bin " << idx << " vs " << bins[idx].pointSize() << " total rP in bin" << endl;
//         // for (int match_idx : angle_matches) if ((idx2distance[match_idx] <= query_distance+dCut)&&(idx2distance[match_idx] >= query_distance-dCut)) result.push_back(match_idx);
//         result.insert(result.end(),angle_matches.begin(),angle_matches.end());
//         cout << result.size() << " distance matches in bin " << idx << endl;
//     }
//     return result;
// }

/* --- --- --- --- --- findResPairs --- --- --- --- --- */

findResPairs::findResPairs(string resPairDBPath, mstreal _dCut, mstreal _angleCut, mstreal _rmsdCut) : DB(resPairDBPath), dCut(_dCut), rmsdCut(_rmsdCut), angleCut(_angleCut) {
    cout << "Loading residue pair backbone atom distances into APV" << endl;
    mstreal maxDistance = 0;
    boundingBox cis_distanceBB(dCut);
    boundingBox cic_distanceBB(dCut);
    while (DB.hasNext()) {
        resPair* rP = DB.next();
        if (rP->getCaDistance() <= 0) {
            delete rP;
            continue;
        }
        int distance_in_chain = rP->getDistanceInChain();
        if (distance_in_chain < 0) {
            // contact
            cis_ResPairs.push_back(rP);
            cis_distanceBB.update(rP->getbbAtomDistances());
        } else {
            if (distance_in_chain > max_distance_in_chain) {
                delete rP;
                continue;
            }
            // close in chain
            cic_ResPairs.push_back(rP);
            cic_distanceBB.update(rP->getbbAtomDistances());
        }
    }
    cout << "There are " << cis_ResPairs.size() << " contact res pairs and " << cic_ResPairs.size() << " close-in-chain res pairs in the DB total" << endl;
    // resPairMap = new distance3AngleTable(maxDistance,dCut,angleCut);
    // for (int i = 0; i < allResPairs.size(); i++) {
    //     resPair* rP = allResPairs[i];
    //     resPairMap->addResPairToTable(rP,i);
    // }
    int cis_Nbuckets = int(ceil(max(max((cis_distanceBB.getXWidth()), (cis_distanceBB.getYWidth())), (cis_distanceBB.getZWidth()))/0.5)); //0.5 Å is the characteristic distance, eq. copied from ProximitySearch MST class
    cis_PS = ProximitySearch(cis_distanceBB.getXMin(),cis_distanceBB.getYMin(),cis_distanceBB.getZMin(),cis_distanceBB.getXMax(),cis_distanceBB.getYMax(),cis_distanceBB.getZMax(),cis_Nbuckets);
    for (int i = 0; i < cis_ResPairs.size(); i++) cis_PS.addPoint(cis_ResPairs[i]->getbbAtomDistances(),i);
    cout << "Done loading contact residue pair information into hash table." << endl; 

    int cic_Nbuckets = int(ceil(max(max((cic_distanceBB.getXWidth()), (cic_distanceBB.getYWidth())), (cic_distanceBB.getZWidth()))/0.5)); //0.5 Å is the characteristic distance, eq. copied from ProximitySearch MST class
    cic_PS = ProximitySearch(cic_distanceBB.getXMin(),cic_distanceBB.getYMin(),cic_distanceBB.getZMin(),cic_distanceBB.getXMax(),cic_distanceBB.getYMax(),cic_distanceBB.getZMax(),cic_Nbuckets);
    for (int i = 0; i < cic_ResPairs.size(); i++) cic_PS.addPoint(cic_ResPairs[i]->getbbAtomDistances(),i);
    cout << "Done loading close-in-chain residue pair information into hash table." << endl; 
}

int findResPairs::searchForMatches(bool verbose) {
    // step 1: quickly eliminate residue pairs by looking at backbone atom distances
    cout << "Backbone atom distance cutoff used for searching: " << dCut << " Å" << endl;
    // cout << "Angle cutoff used for searching: " << angleCut << " degrees" << endl;

    // reset from previous search
    matches.clear();
    verifiedMatches.clear();

    mstreal queryDistance = queryRP.getCaDistance();
    CartesianPoint queryBBdistances = queryRP.getbbAtomDistances();
    // CartesianPoint queryAngles = queryRP.getAngles();
    // cout << "Ca distance: " << queryDistance << endl;
    // cout << "x axis angle: " << queryAngles.getX() << endl;
    // cout << "y axis angle: " << queryAngles.getY() << endl;
    // cout << "z axis angle: " << queryAngles.getZ() << endl;

    // search for residue pairs with matching distances
    timer.start();
    if ((query_distance_in_chain != -1)&&(query_distance_in_chain <= max_distance_in_chain)) {
        cout << "Searching close-in-chain residue pair..." << endl;
        vector<int> matchingResPairIDs = cic_PS.getPointsWithin(queryBBdistances,0,dCut,true);
        // vector<int> matchingResPairIDs = resPairMap->searchResPair(&queryRP);
        int wrong_chain_distance = 0;
        for (int ID : matchingResPairIDs) {
            resPair* rP = cic_ResPairs[ID];
            if (rP->getDistanceInChain() != query_distance_in_chain) {
                wrong_chain_distance++;
                continue;
            }
            matches.push_back(rP);
        }
        timer.stop();
        cout << "Found " << matches.size() << " matching residue pairs by comparing Ca distances in " << timer.getDuration() << " s" << endl;
        cout << "Ignored " << wrong_chain_distance << " potential matches based on distance in chain" << endl;
    } else {
        cout << "Searching contact residue pair..." << endl;
        vector<int> matchingResPairIDs = cis_PS.getPointsWithin(queryBBdistances,0,dCut,true);
        // vector<int> matchingResPairIDs = resPairMap->searchResPair(&queryRP);
        int wrong_chain_distance = 0;
        for (int ID : matchingResPairIDs) {
            resPair* rP = cis_ResPairs[ID];
            matches.push_back(rP);
        }
        timer.stop();
        cout << "Found " << matches.size() << " matching residue pairs by comparing Ca distances in " << timer.getDuration() << " s" << endl;
    }

    // step 2: optimally superimpose each candidate match to the query and calculate RMSD
    if (matches.empty()) return 0;
    vector<Atom*> queryAtoms = queryRP.getAtoms();
    timer.start();
    for (resPair* rP : matches) {
        vector<Atom*> bbAtoms = rP->getAtoms();
        mstreal RMSDval = calc.bestRMSD(bbAtoms,queryAtoms);
        // cout << d << endl;
        // cout << "RMSD: " << RMSDval << endl;
        if (RMSDval <= rmsdCut) {
            // cout << "RMSD: " << RMSDval << endl;
            verifiedMatches.push_back(rP);
        }
    }
    timer.stop();
    cout << "Verified " << verifiedMatches.size() << " residue pairs with RMSD <= " << rmsdCut << " to the query in " << timer.getDuration() << " s" << endl;
    return verifiedMatches.size();
}

int findResPairs::getNumMatchesWithResidueType(bool Ri) {
    int nMatchesWithResType = 0;
    for (resPair* rP : verifiedMatches) {
        if (Ri && (queryRP.getResIAAIndex() == rP->getResIAAIndex())) nMatchesWithResType++;
        if (!Ri && (queryRP.getResJAAIndex() == rP->getResJAAIndex())) nMatchesWithResType++;
    }
    return nMatchesWithResType;
}

vector<int> findResPairs::getNumMatchesByAAType(bool Ri) {
    vector<int> nMatchesAAType(20,0);
    for (resPair* rP : verifiedMatches) {
        if (Ri) nMatchesAAType[rP->getResIAAIndex()]++;
        else nMatchesAAType[rP->getResJAAIndex()]++;
    }
    return nMatchesAAType;
}