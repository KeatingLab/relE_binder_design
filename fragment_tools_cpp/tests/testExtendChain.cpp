#include "mstfasst.h"
#include "mstfuser.h"
#include "mstsystem.h"
#include "mstoptions.h"

#include "alignframes.h"
#include "residuecontact.h"
#include "residueframe.h"


int main(int argc, char *argv[]) {
    MstOptions op;
    op.setTitle("Extends either the N- or C-terminus of a provided fragment. NOTE: very rough beta-test");
    op.addOption("p", "The structure that will be extended",true);
    op.addOption("fasstDB", "Path to the FASST DB that will be searched to extend the chain");
    op.addOption("terminus", "'N' or 'C' for N- or C-terminal extension",true);
    op.addOption("extend_stride", "During each iteration, extend the chain by this many residues (default 4)",false);
    op.addOption("extend_length", "The number of residues to extend the structure by (default: 12)",false);
    op.setOptions(argc,argv);

    MstTimer timer;
    MstUtils::seedRandEngine(42);
    
    Structure S(op.getString("p"));
    string structure_name = MstSys::splitPath(S.getName(),1);
    string dbPath = op.getString("fasstDB");
    string terminus = op.getString("terminus");
    int extend_stride = op.getInt("extend_stride",4);
    int extend_length = op.getInt("extend_length",12);

    cout << "terminus: " << terminus << endl;
    if ((terminus != "C") && (terminus != "N")) MstUtils::error("Must specify either 'N' or 'C' for --terminus");

    FASST F;
    timer.start();
    cout << "Loading fasstdb..." << endl;
    F.readDatabase(dbPath,1);
    timer.stop();
    cout << "Took " << timer.getDuration() << "s to load the DB" << endl;
    fasstSearchOptions options;
    options.setRMSDCutoff(0.3);
    options.setMinNumMatches(1);
    options.setMaxNumMatches(1000);

    fusionParams fParams;
    fusionOutput fOutput;

    fstream pdbOut;
    MstUtils::openFile(pdbOut,"extension.traj.pdb",fstream::out);

    /*
    1) Search the terminus of the structure
    2) Select a match at random and use it to extend to chain
    3) Fuse the matched fragment + extension and go back to 1)
    */
    int starting_structure_length = S.residueSize();
    int current_structure_length = S.residueSize();
    int count = 1;
    while (current_structure_length - starting_structure_length < extend_length) {
        // search
        int resIdxStart, resIdxEnd;
        if (terminus == "N") {
            resIdxStart = 0;
            resIdxEnd = extend_stride;
        } else {
            resIdxStart = current_structure_length - extend_stride;
            resIdxEnd = current_structure_length;
        }
        vector<Residue*> structure_residues = S.getResidues();
        vector<Residue*> query_residues(structure_residues.begin()+resIdxStart,structure_residues.begin()+resIdxEnd);
        Structure query(query_residues);
        F.setQuery(query);
        timer.start();
        fasstSolutionSet sols = F.search();
        timer.stop();
        cout << "Took " << timer.getDuration() << "s to search the DB. Found " << sols.size() << " matches" << endl;
        
        //select
        // clashChecker checker(S);

        Structure extensionFragment;
        Structure extensionOnlyFragment;
        sols.orderByDiscovery(); //not actually random, but good enough for this test
        for (int i = 0; i < sols.size(); i++) {
            fasstSolution selectedSol = sols[i];
            Structure wholeMatch = F.getMatchStructure(selectedSol,false,FASST::matchType::FULL,true);
            vector<int> matchResIdx = F.getMatchResidueIndices(selectedSol);
            Residue* alignmentRes = &wholeMatch.getResidue(matchResIdx.front());
            int resIdxInChain = alignmentRes->getResidueIndexInChain();
            Chain* C = alignmentRes->getChain();
            vector<Residue*> matchChainResidues = C->getResidues();
            int startResIdx, endResIdx;
            int startResIdxExt, endResIdxExt;
            if (terminus == "N") {
                if (resIdxInChain - extend_stride < 0) continue;
                startResIdx = resIdxInChain - extend_stride;
                endResIdx = resIdxInChain + extend_stride + 1;
                startResIdxExt = resIdxInChain - extend_stride + 1;
                endResIdxExt = resIdxInChain - 1;
            } else {
                if (resIdxInChain + 2*extend_stride > C->residueSize() - 1) continue;
                startResIdx = resIdxInChain;
                endResIdx = resIdxInChain + 2*extend_stride;
                startResIdxExt = resIdxInChain + extend_stride + 2;
                endResIdxExt = resIdxInChain + 2*extend_stride - 1;
                cout << startResIdx << " " << endResIdx << endl;
            }
            vector<Residue*> selectedExtensionResidues(matchChainResidues.begin()+startResIdx,matchChainResidues.begin()+endResIdx);
            extensionFragment = Structure(selectedExtensionResidues);
            vector<Residue*> onlyExtensionResidues(matchChainResidues.begin()+startResIdxExt,matchChainResidues.begin()+endResIdxExt);
            // if (checker.checkForClashes(onlyExtensionResidues)) {
            //     cout << "proposed extension clashes with the existing structure" << endl;
            //     continue;
            // }
            break;
        }
        if (extensionFragment.residueSize() == 0) MstUtils::error("Could not find a suitable match to extend from");
        extensionFragment.writePDB("Extension_"+MstUtils::toString(count)+".pdb");

        //fuse
        fusionTopology fTopology(current_structure_length+extend_stride);

        vector<int> baseFragResIdx;
        vector<int> extFragResIdx;
        if (terminus == "N") {
            for (int i = 0; i < current_structure_length; i++) baseFragResIdx.push_back(i+extend_stride);
            for (int i = 0; i < 2*extend_stride; i++) extFragResIdx.push_back(i);
        } else {
            for (int i = 0; i < current_structure_length; i++) baseFragResIdx.push_back(i);
            for (int i = 0; i < 2*extend_stride; i++) extFragResIdx.push_back(i+current_structure_length-extend_stride);     
        }

        fTopology.addFragment(S,baseFragResIdx);
        fTopology.addFragment(extensionFragment,extFragResIdx);

        S = Fuser::fuse(fTopology,fOutput,fParams);
        current_structure_length = S.residueSize();

        // write to file
        pdbOut << "MODEL " << count << endl;
        S.writePDB(pdbOut);
        pdbOut << "ENDMDL" << endl;
        count++;
    }

    return 0;
}