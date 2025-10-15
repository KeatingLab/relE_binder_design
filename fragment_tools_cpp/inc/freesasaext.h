#ifndef _FREESASAEXT_H
#define _FREESASAEXT_H

#include "mstrotlib.h"
#include "msttypes.h"

#include "freesasa.h"

class sasaCalculator {
public:
    sasaCalculator(const Structure& _parentStructure) : parentStructure(_parentStructure) {
        sasaDefined = false;
        params.n_threads = 1;
        setAllowedAA();
        
        prepareSASAStructure();
    };
    
    // sasaCalculator() {
    //     sasaDefined = false;
    //     params.n_threads = 1;
    //     setAllowedAA();
    // };

    ~sasaCalculator() {
        freesasa_structure_free(structure);
        freesasa_node_free(root_node);
    };

    void computeSASA();

//    mstreal getStructureSASA();
    map<Residue*,mstreal> getResidueSASA(bool relative = false);

    map<Atom*,mstreal> getAtomSASA();

protected:
    void setAllowedAA();
    
    void prepareSASAStructure();

    bool sasaDefined = false; //true if sasa has been calculated for current structure
    
    void leftPadTo(string& str, const size_t num, const char paddingChar = ' ');
private:
    //mst variables
    const Structure& parentStructure;
    map<pair<string,int>,Residue*> resMap; //key is chainID and res number
    
    set<string> allowedAA;

    //freesasa variables
    freesasa_parameters params = freesasa_default_parameters;
    freesasa_structure *structure = NULL;
    const freesasa_classifier *classifier = &freesasa_naccess_classifier;
    freesasa_node *root_node = NULL;
};

#endif