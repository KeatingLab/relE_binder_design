#include "searchrespairs.h"

/* --- --- --- --- --- findResPairs --- --- --- --- --- */

findResPairs::findResPairs(string resPairDBPath, mstreal _maxDistance, mstreal _maxRMSD) : DB(resPairDBPath), maxDistance(_maxDistance), maxRMSD(_maxRMSD) {
    cout << "Loading residue pair backbone atom distances into APV" << endl;
    boundingBox distanceBB(1.0);
    while (DB.hasNext()) {
        resPair* rP = DB.next();
        allResPairs.push_back(rP);
        distanceBB.update(rP->getDistances());
    }
    cout << "There are " << allResPairs.size() << " in the DB total" << endl;
    int Nbuckets = int(ceil(max(max((distanceBB.getXWidth()), (distanceBB.getYWidth())), (distanceBB.getZWidth()))/0.5)); //0.5 Å is the characteristic distance, eq. copied from ProximitySearch MST class
    PS = ProximitySearch(distanceBB.getXMin(),distanceBB.getYMin(),distanceBB.getZMin(),distanceBB.getXMax(),distanceBB.getYMax(),distanceBB.getZMax(),Nbuckets);
    for (int i = 0; i < allResPairs.size(); i++) PS.addPoint(allResPairs[i]->getDistances(),i);
    cout << "Done loading into APV." << endl; 
}

int findResPairs::searchForMatches(bool verbose) {
    // step 1: quickly eliminate residue pairs by looking at backbone atom distances
    cout << "distanceCutoff used for searching: " << maxDistance << " Å" << endl;

    // reset from previous search
    matches.clear();
    verifiedMatches.clear();

    CartesianPoint queryDistances = queryRP.getDistances();
    cout << "Ri N - Rj N distance: " << queryDistances.getX() << endl;
    cout << "Ri CA - Rj Ca distance: " << queryDistances.getY() << endl;
    cout << "Ri C - Rj C distance: " << queryDistances.getZ() << endl;

    // search for residue pairs with matching distances
    vector<int> matchingResPairIDs = PS.getPointsWithin(queryDistances,0,maxDistance,true);
    for (int ID : matchingResPairIDs) {
        resPair* rP = allResPairs[ID];
        matches.push_back(rP);
    }
    cout << "Found " << matches.size() << " matching residue pairs by comparing backbone atom distances" << endl;

    // step 2: optimally superimpose each candidate match to the query and calculate RMSD
    if (matches.empty()) return 0;
    vector<Atom*> queryAtoms = queryRP.getAtoms();
    for (resPair* rP : matches) {
        vector<Atom*> bbAtoms = rP->getAtoms();
        mstreal RMSDval = calc.bestRMSD(bbAtoms,queryAtoms);
        // cout << d << endl;
        // cout << "RMSD: " << RMSDval << endl;
        if (RMSDval <= maxRMSD) {
            // cout << "RMSD: " << RMSDval << endl;
            verifiedMatches.push_back(rP);
        }
    }
    cout << "Verified " << verifiedMatches.size() << " residue pairs with RMSD <= " << maxRMSD << " to the query" << endl;
    return verifiedMatches.size();
}

int findResPairs::getNumMatchesWithResidueType(bool Ri) {
    int nMatchesWithResType = 0;
    for (resPair* rP : verifiedMatches) {
        if (Ri && (queryRP.getResIAAIndex() == rP->getResIAAIndex())) nMatchesWithResType++;
        else if (queryRP.getResJAAIndex() == rP->getResJAAIndex()) nMatchesWithResType++;
    }
    return nMatchesWithResType;
}