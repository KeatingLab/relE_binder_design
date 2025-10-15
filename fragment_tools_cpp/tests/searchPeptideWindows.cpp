#include "mstfasst.h"
#include "mstoptions.h"
#include "mstrotlib.h"
#include "mstsystem.h"
#include "msttypes.h"

#include "utilities.h"

int main (int argc, char *argv[]) {
    MstOptions op;
    op.setTitle("Defines overlapping windows of a peptide and searches for matches in a FASST database");
    op.addOption("complexPDB","The path to peptide-protein complex PDB file",true);
    op.addOption("peptideChainID","The peptie chain ID",true);
    op.addOption("fasstDB","The path to a FASST database",true);
    op.addOption("RMSDCutoff","The RMSD cutoff used to define a match (default 0.5)",false);
    op.addOption("segmentLength","The residue length of the overlapping windows (default 6)",false);
    op.setOptions(argc,argv);

    string complexPDBPath = op.getString("complexPDB");
    string peptideChainID = op.getString("peptideChainID");
    string fasstDBPath = op.getString("fasstDB");
    mstreal RMSDCutoff = op.getReal("RMSDCutoff",0.5);
    int segmentLength = op.getInt("segmentLength",6);

    Structure complex(complexPDBPath);
    Chain* peptideC = complex.getChainByID(peptideChainID);
    if (peptideC == NULL) MstUtils::error("Peptide chain not found","searchPeptideWindows::main");
    vector<Atom*> peptideBBAtoms = MiscTools::getBackboneAtoms(peptideC);

    FASST F;
    fasstSearchOptions fOps = F.options();
    cout << "Reading fasst DB..." << endl;
    F.readDatabase(fasstDBPath,1);
    cout << "Done reading DB" << endl;

    string complexName = MstSys::splitPath(complexPDBPath,1);
    string fasstDBName = MstSys::splitPath(fasstDBPath,1);

    fstream data_out;
    MstUtils::openFile(data_out,"peptideWindows.csv",fstream::out);
    data_out << "complex,fasstDB,peptidePosition,segmentLength,numMatches" << endl;

    for (int i = 0; i < peptideC->residueSize() - segmentLength + 1; i++) {
        // Define window of peptide and search against DB
        vector<Atom*> peptideWindow(peptideBBAtoms.begin()+(i*4),peptideBBAtoms.begin()+((i+segmentLength)*4));
        Structure query(peptideWindow);
        F.setQuery(query);
        fOps.setRMSDCutoff(RMSDCutoff);
        F.search();
        
        // Write info to file
        data_out << complexName << "," << fasstDBName << ",";
        data_out << i << "," << segmentLength << "," << F.numMatches();
        data_out << endl;
    }

    cout << "Done!" << endl;
    return 0;
}