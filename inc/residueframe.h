#ifndef _RESIDUEFRAME_H
#define _RESIDUEFRAME_H

#include "mstsequence.h"
#include "msttypes.h"
#include "msttransforms.h"

class residueFrame : public Frame {
    public:
        residueFrame() {};
        residueFrame(Residue* R) : parent(R), aa(SeqTools::aaToIdx(R->getName())) {defineFrame(R);}

        residueFrame* frameRelativeToOther(const residueFrame& other);
        
        CartesianPoint getXPos() {return getO()+getX();}
        CartesianPoint getYPos() {return getO()+getY();}
        CartesianPoint getZPos() {return getO()+getZ();}

        Residue* getParent() {return parent;}
        res_t getAAIdx() {return aa;}

        void setAAIdx(res_t _aa) {aa = _aa;}
        void setReferenceResidue(Residue* reference);

        void writeToFile(string name, fstream& out);

    protected:
        void defineFrame(Residue* R);

    private:
        Residue* parent = nullptr;
        res_t aa = SeqTools::unknownIdx();
};

class augmentedStructure : public Structure {
    public:
        augmentedStructure(string structure_path) : Structure(structure_path,"SKIPHETERO") {defineFrames();}
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