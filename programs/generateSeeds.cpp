#include "msttypes.h"
#include "mstoptions.h"

#include "generateseeds.h"
#include "residuecontact.h"
#include "utilities.h"


int main(int argc, char *argv[]) {
    MstOptions op;
    op.setTitle("Given a target, generate protein backbone fragments around the surface");
    op.addOption("targetPDB","A PDB file containing the structure of the protein target. Must have a defined sequence.",false);
    op.addOption("complexPDB","A PDB file containg structure of a complex between the protein target and some binder. The target must have a defined sequence",false);
    op.addOption("binderChains","--complexPDB mode only. The binder chain ID(s) delimited by '_', e.g. '0_1'",false);
    op.addOption("wholeSurface","--complexPDB mode only. If specified, will generate seeds around the whole surface, not just the site that the binder is interacting with",false);
    op.addOption("numMatches","Find up to this many matches per binding site fragment (if more are found, take the lowest RMSD matches) (default = 1000)",false);
    op.addOption("fasstDB","A path to a structure database that is compatible with FASST",true);
    op.addOption("seqConst","If provided, constrain the matches to the binding site to those with the same central amino acid");
    op.setOptions(argc,argv);

    if (!(op.isGiven("targetPDB"))&&!(op.isGiven("complexPDB") && op.isGiven("binderChains"))) MstUtils::error("Must provide either --targetPDB or --complexPDB and --binderChains");

    string targetPDB = op.getString("targetPDB","");
    string complexPDB = op.getString("complexPDB","");
    string binderChainsString = op.getString("binderChains","");
    string fasstDBPath = op.getString("fasstDB");
    int numMatches = op.getInt("numMatches",1000);
    bool wholeSurface = op.isGiven("wholeSurface");
    bool seqConst = op.isGiven("seqConst");

    seedGenParams params; 
    params.maxNumMatches = numMatches;
    params.seqConst = seqConst;

    // generate seeds
    string binPath, targetName;
    if (targetPDB != "") {
        Structure target = Structure(targetPDB,"SKIPHETERO|ALLOW ILE CD1");

        seedGenerator seedGen(target,fasstDBPath,params);
        binPath = seedGen.generateSeeds();
        targetName = seedGen.getTargetName();
    } else {
        Structure complex(complexPDB,"SKIPHETERO|ALLOW ILE CD1");
        Structure target = complex;

        // get binder/target chains
        vector<string> binderChainIDs = MstUtils::split(binderChainsString);
        vector<Chain*> targetChains, binderChains;
        for (int i = 0; i < target.chainSize(); i++) {
            Chain* C = &target.getChain(i);
            if (find(binderChainIDs.begin(),binderChainIDs.end(),C->getID()) != binderChainIDs.end()) {
                binderChains.push_back(C);
            } else {
                targetChains.push_back(C);
            }
        }

        // define binding site by contacts
        vector<Residue*> bindingSiteRes;
        if (!wholeSurface) {
            vdwContacts vdwC(targetChains,binderChains);
            vector<pair<Residue*,Residue*>> interacting = vdwC.getInteractingResPairs();
            for (auto pair : interacting) {
                Residue* targetR = pair.first;
                if (find(bindingSiteRes.begin(),bindingSiteRes.end(),targetR) == bindingSiteRes.end()) bindingSiteRes.push_back(targetR);
            }
            cout << "Defined " << bindingSiteRes.size() << " binding site residues by vdW contacts" << endl;
        }

        // remove binder chain(s)
        for (Chain* C : binderChains) target.deleteChain(C);

        seedGenerator seedGen(target,fasstDBPath,params);

        if (!wholeSurface) {
            seedGen.setBindingSite(bindingSiteRes);
        }

        // generate seeds
        binPath = seedGen.generateSeeds();
        targetName = seedGen.getTargetName();
        seedGen.writeBindingSiteFragments(targetName+"_bindingSiteFragments.pdb");

        // get the RMSD of each seed to the lowest RMSD window of the peptide
        seedBinaryFile seedBin(binPath);
        fstream out;
        MstUtils::openFile(out,targetName+"_peptideRMSD.csv",fstream::out);
        out << "seedName,rmsd,peptideChainID,peptideResIdx,seedResIdx,length" << endl;

        vector<Chain*> binderChainsInComplex;
        for (int i = 0; i < complex.chainSize(); i++) {
            Chain* C = &complex.getChain(i);
            if (find(binderChainIDs.begin(),binderChainIDs.end(),C->getID()) != binderChainIDs.end()) {
                binderChainsInComplex.push_back(C);
            }
        }
        cout << "Found " << binderChainsInComplex.size() << " binder chains" << endl;
        while (seedBin.hasNext()) {
            Structure* seed = seedBin.next();
            Chain* seedChain = seed->getChainByID("0");
            if (seedChain == NULL) MstUtils::error("Seed missing chain 0","generateSeeds::main");

            MiscTools::alignment bestAlignment;
            string bestAlignmentChainID = "";
            for (Chain* binderChain : binderChainsInComplex) {
                MiscTools::alignment newAlignment = MiscTools::bestRMSD(binderChain,seedChain);
                if (newAlignment.rmsd < bestAlignment.rmsd) {
                    bestAlignment = newAlignment;
                    bestAlignmentChainID = binderChain->getID();
                }
            }
            out << seed->getName() << "," << bestAlignment.rmsd << "," << bestAlignmentChainID << ",";
            out << bestAlignment.CiResIdx << "," << bestAlignment.CjResIdx << "," << bestAlignment.length;
            out << endl;

            delete seed;
        }
        out.close();
    }

    // sample seeds from binary file for visualization in Pymol
    seedBinaryFile seedBin(binPath);
    int numToSample = 500;
    int numSeeds = seedBin.structureCount();
    mstreal sampleRate = min(1.0,mstreal(numToSample)/mstreal(numSeeds));
    cout << "sample rate: " << sampleRate << endl;

    fstream sout;
    MstUtils::openFile(sout,targetName+"_sampledSeeds.pdb",fstream::out);

    while (seedBin.hasNext()) {
        if (MstUtils::randUnit() > sampleRate) {
            seedBin.skip();
            continue;
        }
        Structure* seed = seedBin.next();
        sout << "HEADER    " << seed->getName() << endl;
        seed->writePDB(sout);
    }
    sout.close();

    cout << "Done" << endl;
    return 0;
}