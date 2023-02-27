#include "searchrespairs.h"

/* --- --- --- --- --- distance3AngleTable --- --- --- --- --- */

vector<int> distance3AngleTable::searchResPair(resPair* rP) {
    mstreal query_distance = rP->getCaDistance();
    vector<int> binIdxs = getBinIdxRange(query_distance,dCut);
    vector<int> result;
    for (int idx : binIdxs) {
        vector<int> angle_matches = bins[idx].getPointsWithin(rP->getAngles(),0,angleCut,true);
        cout << angle_matches.size() << " angle matches in bin " << idx << " vs " << bins[idx].pointSize() << " total rP in bin" << endl;
        // for (int match_idx : angle_matches) if ((idx2distance[match_idx] <= query_distance+dCut)&&(idx2distance[match_idx] >= query_distance-dCut)) result.push_back(match_idx);
        result.insert(result.end(),angle_matches.begin(),angle_matches.end());
        cout << result.size() << " distance matches in bin " << idx << endl;
    }
    return result;
}

/* --- --- --- --- --- findResPairs --- --- --- --- --- */

findResPairs::findResPairs(string resPairDBPath, mstreal _dCut, mstreal _angleCut, mstreal _rmsdCut) : DB(resPairDBPath), dCut(_dCut), rmsdCut(_rmsdCut), angleCut(_angleCut) {
    cout << "Loading residue pair backbone atom distances into APV" << endl;
    mstreal maxDistance = 0;
    boundingBox distanceBB(dCut);
    while (DB.hasNext()) {
        resPair* rP = DB.next();
        if (rP->getCaDistance() <= 0) {
            delete rP;
            continue;
        } 
        allResPairs.push_back(rP);
        distanceBB.update(rP->getbbAtomDistances());
    }
    cout << "There are " << allResPairs.size() << " in the DB total" << endl;
    // resPairMap = new distance3AngleTable(maxDistance,dCut,angleCut);
    // for (int i = 0; i < allResPairs.size(); i++) {
    //     resPair* rP = allResPairs[i];
    //     resPairMap->addResPairToTable(rP,i);
    // }
    int Nbuckets = int(ceil(max(max((distanceBB.getXWidth()), (distanceBB.getYWidth())), (distanceBB.getZWidth()))/0.5)); //0.5 Å is the characteristic distance, eq. copied from ProximitySearch MST class
    PS = ProximitySearch(distanceBB.getXMin(),distanceBB.getYMin(),distanceBB.getZMin(),distanceBB.getXMax(),distanceBB.getYMax(),distanceBB.getZMax(),Nbuckets);
    for (int i = 0; i < allResPairs.size(); i++) PS.addPoint(allResPairs[i]->getbbAtomDistances(),i);
    cout << "Done loading residue pair information into hash table." << endl; 
}

int findResPairs::searchForMatches(bool verbose) {
    // step 1: quickly eliminate residue pairs by looking at backbone atom distances
    cout << "Ca distance cutoff used for searching: " << dCut << " Å" << endl;
    cout << "Angle cutoff used for searching: " << angleCut << " degrees" << endl;

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
    vector<int> matchingResPairIDs = PS.getPointsWithin(queryBBdistances,0,dCut,true);
    // vector<int> matchingResPairIDs = resPairMap->searchResPair(&queryRP);
    for (int ID : matchingResPairIDs) {
        resPair* rP = allResPairs[ID];
        matches.push_back(rP);
    }
    timer.stop();
    cout << "Found " << matches.size() << " matching residue pairs by comparing Ca distances and angles in " << timer.getDuration() << " s" << endl;

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