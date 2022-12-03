#include "mstfasst.h"
#include "mstoptions.h"
#include "mstsystem.h"
#include "msttypes.h"

#include "alignframes.h"
#include "scoreinterface.h"
#include "residuecontact.h"
#include "residueframe.h"

int main(int argc, char *argv[]) {
    MstOptions op;
    op.setTitle("Defines protein-protein contacts and finds the prevalence of the pairwise interactions in the PDB");
    op.addOption("structDB","Path to a FASST database of protein structures",true);
    op.addOption("pdb","A PDB file defining the complex between the protein target and designed binder. The target must have a defined sequence",true);
    op.addOption("target","The protein chain ID(s) delimited by '_', e.g. 'A_B'",true);
    op.addOption("binder","The binder chain ID(s) delimited by '_', e.g. '0_1'",true);
    op.addOption("RMSD","The RMSD cutoff (Å) used to define structural matches (default: 0.5)",false);
    // op.addOption("homThresh","The sequence homology threshold that will be applied when determining if a match is redundant. When the value is negative, no homology check is performed (default = 0.6)",false);
    op.setOptions(argc,argv);

    string fasstDBPath = op.getString("structDB");
    string pdbPath = op.getString("pdb");
    string targetChainIDsString = op.getString("target");
    string binderChainIDsString = op.getString("binder");
    mstreal RMSD = op.getReal("RMSD",0.5);

    // get chain ID lists
    vector<string> targetChainIDs = MstUtils::split(targetChainIDsString,"_");
    cout << "target chain IDs: ";
    for (string s : targetChainIDs) cout << " " << s;
    cout << endl;

    vector<string> binderChainIDs = MstUtils::split(binderChainIDsString,"_");
    cout << "binder chain IDs:";
    for (string s : binderChainIDs) cout << " " << s;
    cout << endl;

    // Load the structure and define interface contacts
    string complexName = MstSys::splitPath(pdbPath,1);
    Structure complex(pdbPath);
    // vector<Atom*> complexBBAtoms = MiscTools::getBackboneAtoms(complex);
    // Structure complexBB(complexBBAtoms);

    vector<Chain*> targetChains;
    vector<Chain*> binderChains;
    for (int i = 0; i < complex.chainSize() ; i++) {
        Chain* C = &complex.getChain(i);
        for (string targetChainID : targetChainIDs) {
            if (C->getID() == targetChainID) {
                targetChains.push_back(C);
                continue;
            }
        }
        for (string binderChainID : binderChainIDs) {
            if (C->getID() == binderChainID) {
                binderChains.push_back(C);
                continue;
            }
        }
    }

    vdwContacts conts(targetChains,binderChains);
    vector<pair<Residue*,Residue*>> interactingResPairs = conts.getInteractingResPairs();
    cout << interactingResPairs.size() << " interface contacts when using vdW definition" << endl;

    // Load the FASST DB
    FASST F;
    F.setRMSDCutoff(RMSD);
    cout << "Set FASST RMSD cutoff to: " << F.getRMSDCutoff() << endl; 
    cout << "Reading file as FASST DB" << endl;
    F.readDatabase(fasstDBPath);
    cout << "Done reading DB" << endl;

    // Search each interface contact for matches and report statistics
    fstream info_out;
    MstUtils::openFile(info_out,complexName+"_interfaceCoverage.csv",fstream::out);
    info_out << "complexName,targetRes,targetResName,binderRes,nMatches,nNativeAA" << endl;

    fstream pdb_out;
    MstUtils::openFile(pdb_out,complexName+"_interfaceResiduePairMatches.pdb",fstream::out);

    MstTimer timer;
    for (pair<Residue*,Residue*> resPair : interactingResPairs) {
        // set query
        Structure query(vector<Residue*>({resPair.first,resPair.second}));
        F.setQuery(query);

        // search for matches
        cout << "Searching " << resPair.first->getChainID() << resPair.first->getNum()  << " and " << resPair.second->getChainID() << resPair.second->getNum() << endl;
        timer.start();
        fasstSolutionSet matches = F.search();
        timer.stop();
        cout << "Took " << timer.getDuration() << " s to find " << matches.size() << " matches"<< endl;

        info_out << complexName << ",";
        info_out << resPair.first->getChainID() << resPair.first->getNum() <<  "," << resPair.first->getName() << ",";
        info_out << resPair.second->getChainID() << resPair.second->getNum() << ",";
        info_out << matches.size() << ",";

        int nNativeAA = 0;
        string nativeAAname = resPair.first->getName();
        for (fasstSolution sol : matches) {
            Structure match = F.getMatchStructure(sol,true);
            if (match.getResidue(1).getName() == nativeAAname) nNativeAA++;
        }
        info_out << nNativeAA << endl;

        if (matches.size() > 0) {
            int N = min(5,int(matches.size()));
            for (int i = 0; i < N; i++) {
                fasstSolution sol = matches[i];
                Structure match = F.getMatchStructure(sol,true);
                pdb_out << "HEADER     " << resPair.first->getChainID() << resPair.first->getNum()  << "_" << resPair.second->getChainID() << resPair.second->getNum() << i << endl;
                match.writePDB(pdb_out);
            }
        }
    }

    cout << "Done" << endl;
    return 0;
}