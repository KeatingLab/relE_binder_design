#include "hashframes.h"

/* --- --- --- --- --- boundingBox --- --- --- --- --- */

void boundingBox::update(vector<mobileFrame*> frames) {
    for (mobileFrame* frame : frames) update(frame);
}

void boundingBox::update(Frame* frame) {
    update(frame->getO());
}

bool boundingBox::isWithinBounds(Frame* frame, mstreal tolerance) const {
    return isWithinBounds(frame->getO());
}

bool boundingBox::isWithinBounds(Atom* A, mstreal tolerance) const {
    return isWithinBounds(A->getCoor());
}

bool boundingBox::isWithinBounds(const CartesianPoint& point, mstreal tolerance) const {
    const mstreal& x = point[0];
    const mstreal& y = point[1];
    const mstreal& z = point[2];

    if (x < xMin - tolerance) return false;
    if (x > xMax + tolerance) return false;
    if (y < yMin - tolerance) return false;
    if (y > yMax + tolerance) return false;
    if (z < zMin - tolerance) return false;
    if (z > zMax + tolerance) return false;
    return true;
}

void boundingBox::update(const CartesianPoint& point) {
    const mstreal& x = point[0];
    const mstreal& y = point[1];
    const mstreal& z = point[2];

    if (x - pad < xMin) xMin = x - pad;
    if (x + pad > xMax) xMax = x + pad;
    if (y - pad < yMin) yMin = y - pad;
    if (y + pad > yMax) yMax = y + pad;
    if (z - pad < zMin) zMin = z - pad;
    if (z + pad > zMax) zMax = z + pad;
}

void boundingBox::printBounds() {
    cout << "New bounding box with dimensions: ";
    cout << " x_min = " << xMin << ", x_max = " << xMax;
    cout << " y_min = " << yMin << ", y_max = " << yMax;
    cout << " z_min = " << zMin << ", z_max = " << zMax;
    cout << endl;
}

/* --- --- --- --- --- positionHasher --- --- --- --- --- */

positionHasher::positionHasher(boundingBox _bbox, mstreal _increment) : bbox(_bbox), increment(_increment) {
    if ((increment <= 0) || (increment > 5)) MstUtils::error("increment should be greater than 0 and less than 5","positionHasher::positionHasher()");
    numXBins = ceil((bbox.getXWidth()) / increment);
    numYBins = ceil((bbox.getYWidth()) / increment);
    numZBins = ceil((bbox.getZWidth()) / increment);
    if (getNumBins() > INT_MAX) MstUtils::error("Too many bins (integer overflow). Max allowable value is "+MstUtils::toString(INT_MAX));
}

int positionHasher::hashFrame(Frame* frame, bool strict) {
    const CartesianPoint& origin = frame->getO();
    return hashPoint(origin);
}

int positionHasher::hashPoint(const CartesianPoint& point, bool strict) {
    if (point.size() != 3) MstUtils::error("Hasher only supports 3D coordinates","positionHasher::hashPoint()");
    if (!bbox.isWithinBounds(point)) {
        if (strict) MstUtils::error("Point falls outside of bounding box: "+MstUtils::toString(point.getX())+","+MstUtils::toString(point.getY())+","+MstUtils::toString(point.getZ()),"positionHasher::hashPoint");
        else return -1;
    }
    const mstreal& x = point[0]; const mstreal& y = point[1]; const mstreal& z = point[2];
    
    // Find the bin position in terms of x, y, z
    int xBin = floor((x - bbox.getXMin()) / increment);
    int yBin = floor((y - bbox.getYMin()) / increment);
    int zBin = floor((z - bbox.getZMin()) / increment);
    return hashCoordinates(xBin,yBin,zBin);
}

vector<int> positionHasher::region(Frame* frame, mstreal cutoff) {
    const CartesianPoint& origin = frame->getO();
    return region(origin,cutoff);
}

vector<int> positionHasher::region(const CartesianPoint& point, mstreal cutoff) {
    if (point.size() != 3) MstUtils::error("Hasher only supports 3D coordinates","positionHasher::hashPoint");
    const mstreal& x = point[0];
    const mstreal& y = point[1];
    const mstreal& z = point[2];
    
    // Compute bin index bounds in each direction
    int minX = max(0,int(floor((x - cutoff - bbox.getXMin()) / increment)));
    int maxX = min(numXBins-1,int(ceil((x + cutoff - bbox.getXMin()) / increment)));

    int minY = max(0,int(floor((y - cutoff - bbox.getYMin()) / increment)));
    int maxY = min(numYBins-1,int(ceil((y + cutoff - bbox.getYMin()) / increment)));

    int minZ = max(0,int(floor((z - cutoff - bbox.getZMin()) / increment)));
    int maxZ = min(numZBins-1,int(ceil((z + cutoff - bbox.getZMin()) / increment)));

    // Add all hash values
    vector<int> result;
    for (int x = minX; x <= maxX; x++) {
        for (int y = minY; y <= maxY; y++) {
            for (int z = minZ; z <= maxZ; z++) {
                result.push_back(hashCoordinates(x,y,z));
            }
        }
    }
    return result;
}

CartesianPoint positionHasher::getCenterCoordinatesOfBin(int hash) {
    int xBin, yBin, zBin;
    getCoordinatesFromHash(hash, xBin, yBin, zBin);
    return CartesianPoint((bbox.getXMin()+xBin*increment+increment/2),(bbox.getYMin()+yBin*increment+increment/2),(bbox.getZMin()+zBin*increment+increment/2));
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
}

void positionHasher::getCoordinatesFromHash(int hash, int& xBin, int& yBin, int& zBin) {
    zBin = hash / (numXBins * numYBins);
    yBin = (hash % numXBins * numYBins) / numXBins;
    xBin = ((hash % numXBins * numYBins) % numXBins);
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
    if (alphaAngle >= 2*M_PI) alphaAngle = 0.0;

    // beta is angle around Y-axis (in the Z,X plane)
    const mstreal& Z_x = Z[0];
    const mstreal& Z_z = Z[2];
    betaAngle = fmod(atan2(Z_x,Z_z) + 2*M_PI, 2*M_PI);
    if (betaAngle >= 2*M_PI) betaAngle = 0.0;

    // gamma is angle around Z-axis (in the X,Y plane)
    const mstreal& X_y = X[1];
    const mstreal& X_x = X[0];
    gammaAngle = fmod(atan2(X_y,X_x) + 2*M_PI, 2*M_PI);
    if (gammaAngle >= 2*M_PI) gammaAngle = 0.0;
}

bool anglesFromFrame::L1Check(const anglesFromFrame& other, mstreal cutoff) {
    // This is similar to the L1-Norm, but not the same because angle distance is not a metric
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
    if ((numBins <= 0)) MstUtils::error("Must have at least one bin","orientationHasher::orientationHasher");
    increment = 2*M_PI / numBins;
    if (getNumBinsTotal() > INT_MAX) MstUtils::error("Too many bins (integer overflow). Max allowable value is "+MstUtils::toString(INT_MAX));
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
    getAngleBins(angles,alphaBin,betaBin,gammaBin,false);
    return hashAngleBins(alphaBin,betaBin,gammaBin);
}

void orientationHasher::getAngleBins(const anglesFromFrame& angles, int& alphaBin, int& betaBin, int& gammaBin, bool strict) {
    // perform division in log space for increased precision (directly dividing doubles was not accurate)
    if (strict) {
        alphaBin = int(floor(angles.getAlpha()/increment));
        betaBin = int(floor(angles.getBeta()/increment));
        gammaBin = int(floor(angles.getGamma()/increment));
    } else {
        alphaBin = int(floor(angles.getAlpha()/increment)) % numBins;
        betaBin = int(floor(angles.getBeta()/increment)) % numBins;
        gammaBin = int(floor(angles.getGamma()/increment)) % numBins;
    }
    if ((alphaBin < 0) || (alphaBin >= numBins) || (betaBin < 0) || (betaBin >= numBins) || (gammaBin < 0) || (gammaBin >= numBins)) {
        std::cout << std::fixed << std::setprecision(std::numeric_limits<long double>::digits10 + 1);
        // cout << "2*M_PI: " << 2*M_PI << endl;
        // cout << "increment: " << increment << endl;

        // cout << "alphaAngle: " << angles.getAlpha() << endl;
        // cout << "2*M_PI - alphaAngle: " << 2*M_PI - angles.getAlpha() << endl;
        // cout << "angles.getAlpha() / increment: " << angles.getAlpha() / increment << endl;
        // cout << "exp(log(angles.getAlpha()) - log(increment)): " << exp(log(angles.getAlpha()) - log(increment)) << endl;
        
        // cout << "betaAngle: " << angles.getBeta() << endl;
        // cout << "2*M_PI - betaAngle: " << 2*M_PI - angles.getBeta() << endl;
        // cout << "angles.getBeta() / increment: " << angles.getBeta() / increment << endl;
        // cout << "exp(log(angles.getBeta()) - log(increment)): " << exp(log(angles.getBeta()) - log(increment)) << endl;
        
        // cout << "gammaAngle: " << angles.getGamma() << endl;
        // cout << "2*M_PI - gammaAngle: " << 2*M_PI - angles.getGamma() << endl;
        // cout << "angles.getGamma() / increment: " << angles.getGamma() / increment << endl;
        // cout << "exp(log(angles.getGamma()) - log(increment)): " << exp(log(angles.getGamma()) - log(increment)) << endl;
        MstUtils::error("Provided angles: "+MstUtils::toString(angles.getAlpha())+","+MstUtils::toString(angles.getBeta())+","+MstUtils::toString(angles.getGamma())+" result in invalid bin values: "+MstUtils::toString(alphaBin)+","+MstUtils::toString(betaBin)+","+MstUtils::toString(gammaBin),"orientationHasher::hashAngleBins()");
    }
}

int orientationHasher::hashAngleBins(const int& alphaBin, const int& betaBin, const int& gammaBin) {
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
    allFrames.push_back(frame);
    int posHash = posHasher.hashFrame(frame);
    int rotHash = rotHasher.hashFrame(frame);
    posLookupTable[posHash].push_back(frame);
    rotLookupTable[rotHash].push_back(frame);
}

int frameTable::countSimilarFrames(Frame* frame, mstreal distCut, mstreal angCut, vector<bool>* posHashToIgnore) {
    return findSimilarFrames(frame,distCut,angCut,posHashToIgnore).size();
}

set<mobileFrame*> frameTable::findSimilarFrames(Frame* frame, mstreal distCut, mstreal angCut, vector<bool>* posHashToIgnore) {
    mstreal radCut = M_PI*angCut/(180.0); //rot hasher uses radians
    cout << "Searching for matches to query frame: ";
    printFrameInfo(frame);
    vector<int> posNeigbourHashes = posHasher.region(frame,distCut);
    vector<int> rotNeighborHashes = rotHasher.region(frame,radCut);
    vector<mobileFrame*> posNeighbors; vector<mobileFrame*> rotNeighbors;
    for (int hash : posNeigbourHashes) {
        if ((posHashToIgnore != nullptr) && (posHashToIgnore->at(hash))) continue;
        vector<mobileFrame*> frames = posLookupTable[hash];
        if (!frames.empty()) posNeighbors.insert(posNeighbors.end(),frames.begin(),frames.end());
    }
    // cout << "Found " << posNeigbourHashes.size() << " frames within distance cutoff" << endl;
    for (int hash : rotNeighborHashes) {
        vector<mobileFrame*> frames = rotLookupTable[hash];
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

/* --- --- --- --- --- frameProbability --- --- --- --- --- */

void frameProbability::setTargetResidue(int resIdx, const Structure& target, bool clashCheck) {
    if (targetResDefined.find(resIdx) != targetResDefined.end()) MstUtils::error("A target residue with this index ("+MstUtils::toString(resIdx)+") has already been set","frameProbability::setTargetResidue");
    Structure targetCopy(target);
    Residue& targetRes = target.getResidue(resIdx);
    occupiedVolMap[resIdx].resize(posHasher.getNumBins(),false);
    if (!clashCheck) return;

    MstTimer timer; timer.start();

    // align targetRes to global reference frame and apply transformation to the whole structure
    Transform targetToGlobal = TransformFactory::switchFrames(Frame(),residueFrame(&targetRes));
    targetToGlobal.apply(targetCopy);

    // get the atoms of the structure within the bounding box (+/- some tolerance)
    const boundingBox& bbox = getBoundingBox();
    vector<Atom*> atomsInBox;

    for (Atom* A : targetCopy.getAtoms()) {
        mstreal maxVDWRad = atomVDWcheck.maxSumRadii();
        if (bbox.isWithinBounds(A,maxVDWRad)) atomsInBox.push_back(A);
    }
    cout << atomsInBox.size() << " target atoms in the bounding box" << endl;

    // for each atom, assign the corresponding position bins as occupied
    mstreal CaRadius = atomVDWcheck.getRadius("ALA","CA");
    for (Atom* A : atomsInBox) {
        CartesianPoint point = A->getCoor();
        mstreal distanceToCheck = clashTol*(atomVDWcheck.getRadius(A) + CaRadius);
        vector<int> neighborhood = posHasher.region(point,distanceToCheck);
        for (int bin : neighborhood) {
            // get bin coordinates, check if within distance
            CartesianPoint binCenter = posHasher.getCenterCoordinatesOfBin(bin);
            if (point.distance(binCenter) < distanceToCheck) {
                occupiedVolMap[resIdx][bin] = true;
            }
        }
    }

    // sum up the counts in the remaining bins
    normalizingConstant[resIdx] = 0;
    int numOccupiedBins = 0;
    for (int bin = 0; bin < posHasher.getNumBins(); bin++) {
        if (!occupiedVolMap[resIdx][bin]) {
            normalizingConstant[resIdx] += getNumFramesInPosBin(bin);
        } else {
            numOccupiedBins++;
        }

    }
    timer.stop();
    cout << "Took " << timer.getDuration() << " s to align structure to target residue and label voxels corresponding to occupied volume" << endl;
    cout << numOccupiedBins << " bins out of " << posHasher.getNumBins() << " total are occupied" << endl;
    cout << "Normalizing constant for target residue with index " << resIdx << " is " << normalizingConstant[resIdx] << " out of " << getTotalNumFramesInTable() << endl;
    targetResDefined.insert(resIdx);
}

pair<int,int> frameProbability::findInteractingFrameProbability(Frame* frame, int targetResIdx) {
    // check if frame overlaps with the volume occupied by the protein
    int hash = posHasher.hashFrame(frame);
    // if (occupiedVolMap[targetResIdx][hash]) return pair<int,int>(-1.0,normalizingConstant[targetResIdx]);

    // find the number of matches to the frame in the DB and normalize
    int numMatches = countSimilarFrames(frame, distCut, oriCut, &(occupiedVolMap[targetResIdx]));
    return pair<int,int>(numMatches,normalizingConstant[targetResIdx]);
}

vector<mstreal> frameProbability::findBinderResAADist(Frame* frame, int targetResIdx, mstreal pseudocount) {
    // check if frame overlaps with the volume occupied by the protein
    int hash = posHasher.hashFrame(frame);

    // find the number of matches to the frame in the DB and normalize
    set<mobileFrame*> matches = findSimilarFrames(frame, distCut, oriCut, &(occupiedVolMap[targetResIdx]));
    return aaDistFromMobileFrames(matches,pseudocount);
}

vector<mstreal> frameProbability::aaDistFromMobileFrames(set<mobileFrame*> frames, mstreal pseudocount) {
    vector<mstreal> aaDist;
    for (string aa : SeqToolsExtension::getAANames()) {
        res_t aaIdx = SeqTools::aaToIdx(aa);
        mstreal count = pseudocount;
        for (mobileFrame* frame : frames) {
            if (frame->getResJIndex() == aaIdx) count++;
        }
        aaDist.push_back(count);
    }
    return aaDist;
}