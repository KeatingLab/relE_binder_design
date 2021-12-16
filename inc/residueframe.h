#ifndef _RESIDUEFRAME_H
#define _RESIDUEFRAME_H

#include "msttypes.h"
#include "msttransforms.h"

struct sphericalCoordinate {
    void constructSphericalCoordinates(CartesianPoint P, bool unit_sphere) {
        constructSphericalCoordinatesFromXYZ(P.getX(),P.getY(),P.getZ(),unit_sphere);
    };

    //  https://en.wikipedia.org/wiki/Spherical_coordinate_system#Coordinate_system_conversions
    //  using the physics (ISO 80000-2:2019 convention)
    void constructSphericalCoordinatesFromXYZ(mstreal x, mstreal y, mstreal z, bool unit_sphere) {
        if (unit_sphere) r = 1.0;
        else r = (CartesianPoint(x,y,z)).norm();
        
        theta = acos((x/r));

        if (x > 0) {
            phi = atan((y/x));
        } else if (x < 0) {
            phi = atan((y/x)) + M_PI;
        } else phi = M_PI / 2.0;
    };

    sphericalCoordinate operator-(const sphericalCoordinate& other) const {
        sphericalCoordinate relSphericalCoordinate;

        // NOTE: assumes a "unit sphere" (so the two coordinates differ only by rotation)
        if ((this->r != 1.0)||(other.r != 1.0)) MstUtils::error("Either this.r or other.r are not properly normalized (i.e. != 1): "+MstUtils::toString(this->r)+","+MstUtils::toString(other.r),"sphericalCoordinate:operator-");
        relSphericalCoordinate.r = 1.0;

        // find the change in angle, relative to self
        relSphericalCoordinate.theta = other.theta - this->theta;
        relSphericalCoordinate.phi = other.phi - this->phi;
        
        return relSphericalCoordinate;
    };

    mstreal r = 0;
    mstreal theta = 0;
    mstreal phi = 0;
};

class residueFrame : public Frame {
    public:
        residueFrame() {};
        residueFrame(Residue* R);
        residueFrame(const CartesianPoint& _position, const sphericalCoordinate& _orientation) : position(_position),orientation(_orientation) {return;}

        residueFrame frameRelativeToOther(residueFrame* other);

        const CartesianPoint& getPosition() {return position;}
        const sphericalCoordinate& getOrientation() {return orientation;}

        // CartesianPoint getUPos() {return position+u;}
        // CartesianPoint getBPos() {return position+b;}
        // CartesianPoint getNPos() {return position+n;}

        CartesianPoint getNPos() {return getO()+getX();}
        CartesianPoint getBPos() {return getO()+getY();}
        CartesianPoint getUPos() {return getO()+getZ();}
        

        Residue* getParent() {return parent;}

    protected:
        void defineFrame(Residue* R);

    private:
        CartesianPoint position;
        sphericalCoordinate orientation;

        //vectors storing orientation info
        CartesianPoint u;
        CartesianPoint b;
        CartesianPoint n;

        Residue* parent = nullptr;
};

class augmentedStructure : public Structure {
    public:
        augmentedStructure(string structure_path) : Structure(structure_path) {defineFrames();}
        augmentedStructure(const Structure& S): Structure(S) {
              defineFrames();
        }

        residueFrame* getResidueFrame(int res_idx) {
            if ((res_idx < 0)||(res_idx >= frames.size())) MstUtils::error("Provided value "+MstUtils::toString(res_idx)+" is out of range: (0,"+MstUtils::toString(frames.size()-1),"residueFrame::getResidueFrame");
            return &frames[res_idx];
        }

        void writeToFile(string path_prefix);
    protected:
        void defineFrames() {
            for (Residue* R : getResidues()) {
                frames.emplace_back(residueFrame(R));
            }
        }
    private:
        vector<residueFrame> frames;
};

#endif