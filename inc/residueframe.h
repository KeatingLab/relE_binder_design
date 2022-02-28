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
        augmentedStructure() : Structure() {}

        augmentedStructure(string structure_path, bool verbose = false) : Structure(structure_path,"SKIPHETERO") {
            prepareStructure(verbose);
        }

        augmentedStructure(const Structure& S, bool verbose = false): Structure(S) {
            prepareStructure(verbose);
        }

        residueFrame* getResidueFrame(int res_idx) {
            if ((res_idx < 0)||(res_idx >= frames.size())) MstUtils::error("Provided value "+MstUtils::toString(res_idx)+" is out of range: (0,"+MstUtils::toString(frames.size()-1),"residueFrame::getResidueFrame");
            return &(frames[res_idx]);
        }

        void writeToFile(string path_prefix);
    protected:
        void prepareStructure(bool verbose = false) {
            defineFrames();

            // strip hydrogen atoms
            int count = 0;
            for (Residue* R : getResidues()) {
                for (int i = 0; i < R->atomSize(); i++) {
                    Atom* A = &R->getAtom(i);
                    string atomName = A->getName();
                    if (atomName.empty()) MstUtils::error("Atom from residue "+R->getChainID()+MstUtils::toString(R->getNum())+" has no name","augmentedStructure::prepareStructure");
                    if (atomName.at(0) == 'H') {
                        R->deleteAtom(i);
                        count++;
                    }
                }
            }
            if (verbose) cout << "Deleted " << count << " hydrogen atoms" << endl;
        }

        void defineFrames() {
            for (Residue* R : getResidues()) {
                frames.emplace_back(residueFrame(R));
            }
        }
    private:
        vector<residueFrame> frames;
};

#endif