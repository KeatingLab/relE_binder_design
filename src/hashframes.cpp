#include "hashframes.h"

/* --- --- --- --- --- positionHasher --- --- --- --- --- */

positionHasher::positionHasher(vector<mstreal> _bbox, mstreal _increment) : bbox(_bbox), increment(_increment) {
    MstUtils::assert((bbox.size() == 6),"Bounding box should have six values, but has:"+MstUtils::toString(bbox.size()),"positionHasher::positionHasher");
    numXBins = ceil((bbox[1] - bbox[0]) / increment);
    numYBins = ceil((bbox[3] - bbox[2]) / increment);
    numZBins = ceil((bbox[5] - bbox[4]) / increment);
    if ((numXBins*numYBins*numZBins) > INT_MAX) MstUtils::error("Too many bins (integer overflow). Max allowable value is "+MstUtils::toString(INT_MAX));
}

int positionHasher::hashFrame(Frame* frame) {
    const CartesianPoint& o = frame->getO();
    const mstreal& x = o[0]; const mstreal& y = o[1]; const mstreal& z = o[2];
    
    // Find the bin position in terms of x, y, z
    int xBin = floor((x - bbox[0]) / increment);
    int yBin = floor((y - bbox[2]) / increment);
    int zBin = floor((z - bbox[4]) / increment);
    return hashCoordinates(xBin,yBin,zBin);
}

vector<int> positionHasher::region(Frame* frame, mstreal cutoff) {
    vector<int> result;
    
    const CartesianPoint& o = frame->getO();
    const mstreal& o_x = o[0];
    const mstreal& o_y = o[1];
    const mstreal& o_z = o[2];
    
    // Compute bin index bounds in each direction
    int minX = max(0,int(floor((o_x - cutoff - bbox[0]) / increment)));
    int maxX = min(numXBins-1,int(ceil((o_x + cutoff - bbox[0]) / increment)));

    int minY = max(0,int(floor((o_y - cutoff - bbox[2]) / increment)));
    int maxY = min(numYBins-1,int(ceil((o_y + cutoff - bbox[2]) / increment)));

    int minZ = max(0,int(floor((o_z - cutoff - bbox[4]) / increment)));
    int maxZ = min(numZBins-1,int(ceil((o_z + cutoff - bbox[4]) / increment)));

    // Add all hash values
    for (int x = minX; x <= maxX; x++) {
        for (int y = minY; y <= maxY; y++) {
            for (int z = minZ; z <= maxZ; z++) {
                result.push_back(hashCoordinates(x,y,z));
            }
        }
    }
    return result;
}

void positionHasher::updateBoundingBox(vector<mobileFrame*> frames, vector<mstreal>& bbox, bool verbose) {
    for (Frame* frame : frames) {
        const CartesianPoint& o = frame->getO();
        const mstreal& o_x = o[0];
        const mstreal& o_y = o[1];
        const mstreal& o_z = o[2];

        if (o_x < bbox[0]) bbox[0] = o_x;
        if (o_x > bbox[1]) bbox[1] = o_x;
        if (o_y < bbox[2]) bbox[2] = o_y;
        if (o_y > bbox[3]) bbox[3] = o_y;
        if (o_z < bbox[4]) bbox[4] = o_z;
        if (o_z > bbox[5]) bbox[5] = o_z;
    }
    if (verbose) {
        cout << "New bounding box with dimensions: ";
        cout << "x_min = " << bbox[0] << ", x_max = " << bbox[1];
        cout << " y_min = " << bbox[2] << ", y_max = " << bbox[3];
        cout << " z_min = " << bbox[4] << ", z_max = " << bbox[5];
        cout << endl;
    }
}

int positionHasher::hashCoordinates(const int& xBin, const int& yBin, const int& zBin) {
    /*
     The following operation computes a unique hash, which is essentially an address, for each bin
     in the bounding box. Note that if the bounding box becomes too large relative to increment, it
     is possible that we could run out of unique hashes.

     Consider a tensor with dimensions (x,y,z). Each position in the tensor is indexed by 
     some coordinate, (xBin,yBin,zBin) and corresponds to a "bin" or voxel. To convert the address
     into a single number, we walk through the tensor in a consistent manner: 1) start at the origin
     (0,0,0), 2) walk in the X direction until you hit the last position (numXbins,0,0), then take a step 
     in the Y direction and loop back in the X to (0,1,0), 3) continue this until you reach (numXbins,numYbins,0)
     at which point we take loop back in both X,Y and take a step in the Z direction to (0,0,1), etc.
     */
    return zBin * numXBins * numYBins + yBin * numXBins + xBin;
    // return xBin * _numYBins * _numZBins + yBin * _numZBins + zBin; older method
}

/* --- --- --- --- --- anglesFromFrame --- --- --- --- --- */

anglesFromFrame::anglesFromFrame(Frame* frame) {
    const CartesianPoint& X = frame->getX();
    const CartesianPoint& Y = frame->getY();
    const CartesianPoint& Z = frame->getZ();

    /* Rotation is defined around the x, y, and z axis of the global reference frame
    e.g. a rotation around the x-axis is defined like
    
    z-axis
    ^   . [Y_y,Y_z] of the mobile frame
    |  /
    | / alpha = the angle between the Y-axis of mobile frame and the y-axis of global frame
    |/
     --------- > y-axis (x-axis is pointing out of the screen)

     When Y_z < 0, the angle returned is in the range [-M_PI,0]. The following will convert those angles to the range
     [M_PI,2*M_PI] and not change the angles when Y_z > 0.

    alpha = fmod(atan2(Y_z,Y_y) + 2 * M_PI,2 * MPI)
        
        */
    // alpha is angle around X-axis (in the Y,Z plane)
    const mstreal& Y_z = Y[2];
    const mstreal& Y_y = Y[1];
    alphaAngle = fmod(atan2(Y_z,Y_y) + 2*M_PI, 2*M_PI);

    // beta is angle around Y-axis (in the Z,X plane)
    const mstreal& Z_x = Z[0];
    const mstreal& Z_z = Z[2];
    betaAngle = fmod(atan2(Z_x,Z_z) + 2*M_PI, 2*M_PI);

    // gamma is angle around Z-axis (in the X,Y plane)
    const mstreal& X_y = X[1];
    const mstreal& X_x = X[0];
    gammaAngle = fmod(atan2(X_y,X_x) + 2*M_PI, 2*M_PI);
}

bool anglesFromFrame::L1Check(const anglesFromFrame& other, mstreal cutoff) {
    // This is analogous to the L1-Norm, but not the same because angle distance is not a metric
    mstreal otherAlphaAngle = other.getAlpha();
    mstreal otherBetaAngle = other.getBeta();
    mstreal otherGammaAngle = other.getGamma();

   // If the angle is greater than Pi, we know what we need to measure around the other side of the circle
    mstreal alphaDev = abs(alphaAngle - otherAlphaAngle);
    if (alphaDev > M_PI) alphaDev = abs((alphaAngle > otherAlphaAngle) ? alphaAngle-2*M_PI-otherAlphaAngle : otherAlphaAngle-2*M_PI-alphaAngle);
    mstreal betaDev = abs(betaAngle - otherBetaAngle);
    if (betaDev > M_PI) betaDev = abs((betaAngle > otherBetaAngle) ? betaAngle-2*M_PI-otherBetaAngle : otherBetaAngle-2*M_PI-betaAngle);
    mstreal gammaDev = abs(gammaAngle - otherGammaAngle);
    if (gammaDev > M_PI) gammaDev = abs((gammaAngle > otherGammaAngle) ? gammaAngle-2*M_PI-otherGammaAngle : otherGammaAngle-2*M_PI-gammaAngle);
    
    return ((alphaDev + betaDev + gammaDev) <= cutoff);
}

/* --- --- --- --- --- orientationHasher --- --- --- --- --- */

orientationHasher::orientationHasher(int _numBins) : numBins(_numBins) {
    increment = 2*M_PI / numBins;
    if ((numBins*numBins*numBins) > INT_MAX) MstUtils::error("Too many bins (integer overflow). Max allowable value is "+MstUtils::toString(INT_MAX));
}

int orientationHasher::hashFrame(Frame* frame) {
    anglesFromFrame angles(frame);
    return hashAngles(angles);
}

vector<int> orientationHasher::region(Frame* frame, mstreal cutoff) {
    vector<int> result;
    vector<int> alphaBins;
    vector<int> betaBins;
    vector<int> gammaBins;
    getAngleBinsRange(frame,cutoff,alphaBins,betaBins,gammaBins);

    // Add all hash values
    for (int alphaBin : alphaBins) {
        for (int betaBin : betaBins) {
            for (int gammaBin : gammaBins) {
                result.push_back(hashAngleBins(alphaBin,betaBin,gammaBin));
            }
        }
    }
    return result;
}

int orientationHasher::hashAngles(const anglesFromFrame& angles) {
    int alphaBin = -1;
    int betaBin = -1;
    int gammaBin = -1;
    getAngleBins(angles,alphaBin,betaBin,gammaBin);
    return hashAngleBins(alphaBin,betaBin,gammaBin);
}

void orientationHasher::getAngleBins(const anglesFromFrame& angles, int& alphaBin, int& betaBin, int& gammaBin) {
    alphaBin = floor(angles.getAlpha() / increment);
    betaBin = floor(angles.getBeta() / increment);
    gammaBin = floor(angles.getGamma() / increment);
}

int orientationHasher::hashAngleBins(const int& alphaBin, const int& betaBin, const int& gammaBin) {
    if ((alphaBin < 0) || (alphaBin >= numBins) || (betaBin < 0) || (betaBin >= numBins) || (gammaBin < 0) || (gammaBin >= numBins))
     MstUtils::error("bin values invalid: "+MstUtils::toString(alphaBin)+","+MstUtils::toString(betaBin)+","+MstUtils::toString(gammaBin),"orientationHasher::hashAngleBins");
    /* Same idea as in the position hasher. The three angles still form a 3D space, which is now
    technically a 3-torus, but the fact that there is a discontinuity in the hashes doesn't matter
    as long as they're all unique. */
    return gammaBin * numBins * numBins + betaBin * numBins + alphaBin;
}

void orientationHasher::getAngleBinsRange(const anglesFromFrame& angles, mstreal cutoff, vector<int>& alphaBins, vector<int>& betaBins, vector<int>& gammaBins) {
    /*
    Now we pay for binning the angles into vectors.
    If the angle +/- the cutoff crosses the boundary at 2PI -> 0, we need to loop around the circle.
    */
    mstreal alpha = angles.getAlpha(); 
    mstreal beta = angles.getBeta();
    mstreal gamma = angles.getGamma();
    if (cutoff >= M_PI) {
        // return all bins
        alphaBins.reserve(numBins); betaBins.reserve(numBins); gammaBins.reserve(numBins);
        for (int bin = 0; bin < numBins; bin++) {
            alphaBins.push_back(bin);
            betaBins.push_back(bin);
            gammaBins.push_back(bin);
        }
    } else {
        int alphaBinStart = floor((alpha - cutoff) / increment);
        int alphaBinEnd = floor((alpha + cutoff) / increment);
        alphaBins.reserve(alphaBinEnd-alphaBinStart+1);
        if (alphaBinStart < 0) {
            for (int bin = alphaBinStart + numBins; bin < numBins; bin++) alphaBins.push_back(bin);
            for (int bin = 0; bin <= alphaBinEnd; bin++) alphaBins.push_back(bin);
        } else if (alphaBinEnd >= numBins) {
            // crosses discontinuity, add the segment of the circle before/after discontinuity separately
            for (int bin = alphaBinStart; bin < numBins; bin++) alphaBins.push_back(bin);
            for (int bin = 0; bin <= alphaBinEnd - numBins; bin++) alphaBins.push_back(bin);
        } else {
            for (int bin = alphaBinStart; bin <= alphaBinEnd; bin++) alphaBins.push_back(bin);
        }
        int betaBinStart = floor((beta - cutoff) / increment);
        int betaBinEnd = floor((beta + cutoff) / increment);
        betaBins.reserve(betaBinEnd-betaBinStart+1);
        if (betaBinStart < 0) {
            for (int bin = betaBinStart + numBins; bin < numBins; bin++) betaBins.push_back(bin);
            for (int bin = 0; bin <= betaBinEnd; bin++) betaBins.push_back(bin);
        } else if (betaBinEnd >= numBins) {
            // crosses discontinuity, add the segment of the circle before/after discontinuity separately
            for (int bin = betaBinStart; bin < numBins; bin++) betaBins.push_back(bin);
            for (int bin = 0; bin <= betaBinEnd - numBins; bin++) betaBins.push_back(bin);
        } else {
            for (int bin = betaBinStart; bin <= betaBinEnd; bin++) betaBins.push_back(bin);
        }
        int gammaBinStart = floor((gamma - cutoff) / increment);
        int gammaBinEnd = floor((gamma + cutoff) / increment);
        gammaBins.reserve(gammaBinEnd-gammaBinStart+1);
        if (gammaBinStart < 0) {
            for (int bin = gammaBinStart + numBins; bin < numBins; bin++) gammaBins.push_back(bin);
            for (int bin = 0; bin <= gammaBinEnd; bin++) gammaBins.push_back(bin);
        } else if (gammaBinEnd >= numBins) {
            // crosses discontinuity, add the segment of the circle before/after discontinuity separately
            for (int bin = gammaBinStart; bin < numBins; bin++) gammaBins.push_back(bin);
            for (int bin = 0; bin <= gammaBinEnd - numBins; bin++) gammaBins.push_back(bin);
        } else {
            for (int bin = gammaBinStart; bin <= gammaBinEnd; bin++) gammaBins.push_back(bin);
        }
    }
}

/* --- --- --- --- --- frameTable --- --- --- --- --- */

void frameTable::insertFrame(mobileFrame* frame) {
    int posHash = posHasher.hashFrame(frame);
    int rotHash = rotHasher.hashFrame(frame);
    allFrames.push_back(frame);
    posMap[posHash].push_back(frame);
    rotMap[rotHash].push_back(frame);
}

int frameTable::countSimilarFrames(Frame* frame, mstreal distCut, mstreal angCut) {
    return findSimilarFrames(frame,distCut,angCut).size();
}

set<mobileFrame*> frameTable::findSimilarFrames(Frame* frame, mstreal distCut, mstreal angCut) {
    mstreal radCut = 2*M_PI*angCut/(360.0); //rot hasher uses radians
    cout << "Searching for matches to query frame: ";
    printFrameInfo(frame);
    vector<int> posNeigbourHashes = posHasher.region(frame,distCut);
    vector<int> rotNeighborHashes = rotHasher.region(frame,radCut);
    vector<mobileFrame*> posNeighbors; vector<mobileFrame*> rotNeighbors;
    for (int hash : posNeigbourHashes) {
        vector<mobileFrame*> frames = posMap[hash];
        if (!frames.empty()) posNeighbors.insert(posNeighbors.end(),frames.begin(),frames.end());
    }
    // cout << "Found " << posNeigbourHashes.size() << " frames within distance cutoff" << endl;
    for (int hash : rotNeighborHashes) {
        vector<mobileFrame*> frames = rotMap[hash];
        if (!frames.empty()) rotNeighbors.insert(rotNeighbors.end(),frames.begin(),frames.end());
    }
    // cout << "Found " << rotNeighborHashes.size() << " frames within angle cutoff" << endl;
    return verify(frame,posNeighbors,distCut,rotNeighbors,angCut);
}

set<mobileFrame*> frameTable::verify(Frame* queryFrame, const vector<mobileFrame*>& distNeighbors, mstreal distCut, const vector<mobileFrame*>& rotNeighbors, mstreal angCut) {
    set<mobileFrame*> result;
    set<mobileFrame*> pos; set<mobileFrame*> rot;

    // check that the euclidian distance is within the cutoff
    CartesianPoint qCA = queryFrame->getO();
    for (mobileFrame* frame : distNeighbors) {
        CartesianPoint neighborCA = frame->getO();
        if (qCA.distance(neighborCA) <= distCut) pos.insert(frame);
    }

    // check that the norm of the angle deviations is within the cutoff
    anglesFromFrame qAngles(queryFrame);
    for (mobileFrame* frame : rotNeighbors) {
        anglesFromFrame neighborAngles(frame);
        if (qAngles.L1Check(neighborAngles,angCut)) rot.insert(frame);
    }

    // find the intersection of the two sets
    set_intersection(pos.begin(),pos.end(),rot.begin(),rot.end(),std::inserter(result,result.begin()));
    return result;
}

void frameTable::printFrameInfo(Frame* frame) {
    anglesFromFrame angles(frame);
    cout << "x: " << frame->getO()[0];
    cout << "\ty: " << frame->getO()[1];
    cout << "\tz: " << frame->getO()[2];
    cout << "\talpha: " << 360*angles.getAlpha()/(2*M_PI);
    cout << "\tbeta: " << 360*angles.getBeta()/(2*M_PI);
    cout << "\tgamma: " << 360*angles.getGamma()/(2*M_PI);
    cout << endl;
}