#include "msttypes.h"
#include "mstoptions.h"

#include "alignframes.h"
#include "generateseeds.h"
#include "scoreinterface.h"
#include "residuecontact.h"
#include "residueframe.h"

#include <chrono>

int main(int argc, char *argv[]) {
    MstOptions op;
    op.setTitle("Scores a protein-protein interface by finding the prevalence of pairwise interaction geometries in the PDB");
    op.addOption("frameDB","Path to mobile frame database",true);
    op.addOption("contactData","Path to a JSON file describing probability density of contacts",true);
    op.addOption("complexPDB","A PDB file defining the complex between the protein target and designed binder. The target must have a defined sequence",false);
    op.addOption("binderChains","The binder chain ID(s) delimited by '_', e.g. '0_1'. If in --targetPDB mode, will be removed from the structure before scoring",false);
    op.addOption("targetChains","--complexPDB mode only. The protein chain(s) delimited by '_', e.g. 'A_B'",false);
    op.addOption("targetPDB","A PDB file defining the structure of the protein target. Must have a defined sequence.",false);
    op.addOption("binderList","--targetPDB mode only. A file where each line is the path to a binder structure to be scored",false);
    op.addOption("binderBin","--targetPDB mode only. A seed binary file.",false);
    op.addOption("binderBinSubset","--targetPDB and --binderBin mode only. A list of structure names that will be extracted from the seed binary file and scored",false);
    op.addOption("vdwContacts","If provided, will define the interface using VDW contacts, otherwise will use CB definition. Not compatible with --targetPDB mode",false);
    op.addOption("distanceCutoff","The distance cutoff that is applied when determining whether a putative match has the proper position (default = 0.5 Å)",false);
    op.addOption("orientationCutoff","The orientation cutoff that is applied when determining whether a putative match has the proper orientation (default = 15°)",false);
    op.addOption("dontRenormalizeProbability","By default, potential residue interactions that would clash with atoms of the target backbone are omitted and the conditional probability the interaction is renormalized. If this is provided, will not do this",false);
    op.addOption("verbose","");
    op.setOptions(argc,argv);

    if (op.isGiven("complexPDB")) {
        if (!op.isGiven("targetChains") || !op.isGiven("binderChains")) MstUtils::error("If scoring a complex, must provide '--targetChains' and '--binderChains'");
    } else if (op.isGiven("targetPDB")) {
        if ((!op.isGiven("binderList")) && (!op.isGiven("binderBin"))) MstUtils::error("If providing apo target structure, must also provide '--binderList' or '--binderBin'");
        if (!op.isGiven("binderBin") && op.isGiven("binderBinSubset")) MstUtils::error("Can only extract a subset using --binderBinSubset if --binderBin is provided");
    } else {
        MstUtils::error("Must provide either --complexPDB or --targetPDB");
    }

    string mobileFrameDB = op.getString("frameDB");
    string contactData = op.getString("contactData");
    string complexPDB = op.getString("complexPDB","");
    string targetChainIDsString = op.getString("targetChains");
    string binderChainIDsString = op.getString("binderChains");
    string targetPDB = op.getString("targetPDB");
    string binderListPath = op.getString("binderList","");
    string seedBinPath = op.getString("binderBin","");
    string binderBinSubsetPath = op.getString("binderBinSubset","");
    bool defineVDWContacts = op.isGiven("vdwContacts");
    mstreal posCut = op.getReal("distanceCutoff",0.5);
    mstreal oriCut = op.getReal("orientationCutoff",15.0);
    bool renormalizeProbabilities = !op.isGiven("dontRenormalizeProbability");
    bool verbose = op.isGiven("verbose");

    binderScorer* scorer = nullptr;

    binderScorerParams params;
    params.frameDBPath = mobileFrameDB;
    params.potentialContactsJSONPath = contactData;
    params.posCut = posCut;
    params.oriCut = oriCut;
    params.renormalizeProbabilities = renormalizeProbabilities;
    params.verbose = verbose;

    if (complexPDB != "") {
        augmentedStructure complex(complexPDB);
        scorer = new binderScorer(params,complex,binderChainIDsString,targetChainIDsString);

        if (!defineVDWContacts) {
            scorer->defineTargetBindingSiteResiduesByrSASA();
            scorer->defineInterfaceByPotentialContacts();
        }
        mstreal score = scorer->scoreInterface();
        cout << "Interface score: " << score << endl;
        cout << "Non-designable contacts: " << scorer->countNonDesignableContacts() << endl;

        // write the score counts out to a file
        scorer->writeContactScoresToFile();
        scorer->writeResidueClashesToFile();
        scorer->writeBinderScoresToFile();

        scorer->writeContactPropertyToFile();
        scorer->writeResiduesClashesToSimpleFile();
        scorer->writePSSM();

    } else {
        // score target with multiple designed binders mode
        Structure targetS(targetPDB,"SKIPHETERO|ALLOW ILE CD1");

        if (binderChainIDsString != "") {
            vector<string> binderChainIDs = MstUtils::split(binderChainIDsString,"_");
            for (string chainID : binderChainIDs) {
                Chain* C = targetS.getChainByID(chainID);
                if (C == NULL) MstUtils::error("Chain not found in structure: "+MstUtils::toString(chainID),"binderScorer::main()");
                targetS.deleteChain(C);
            }
            cout << "Deleted " << binderChainIDs.size() << " binder chains from the structure" << endl;
        }
        augmentedStructure target(targetS);
        scorer = new binderScorer(params,target);

        // define the target binding site residues using relSASA
        scorer->defineTargetBindingSiteResiduesByrSASA();

        int count = 1;
        if (binderListPath != "") {
            // load each structure, score, and write results
            vector<string> structuresList = MstUtils::fileToArray(binderListPath);
            cout << "Loaded " << structuresList.size() << " binder structure paths to score with the target" << endl;
            for (string structurePath : structuresList) {
                cout << "Loading " << structurePath << " structure " << count << "/" << structuresList.size() << endl;
                augmentedStructure binder(structurePath);
                scorer->setBinder(&binder);
                scorer->defineInterfaceByPotentialContacts();

                scorer->scoreInterface();

                // write the score counts out to a file
                scorer->writeContactScoresToFile();
                scorer->writeResidueClashesToFile();
                scorer->writeBinderScoresToFile();
                scorer->writeContactPropertyToFile("");
                count++;
            }
        } else {
            seedBinaryFile seedBin(seedBinPath);
            int totalNumSeeds = seedBin.structureCount(); int count = 1;
            cout << "Seed binary file has " << seedBin.structureCount() << " seeds" << endl;
            if (binderBinSubsetPath != "") {
                // make directories for more detailed info files
                string pdbFilePrefix = "seedStructures/";
                string contactFilePrefix = "contactFiles/";
                string pssmFilePrefix = "pssmFiles/";
                MstSys::cmkdir(pdbFilePrefix);
                MstSys::cmkdir(contactFilePrefix);
                MstSys::cmkdir(pssmFilePrefix);

                vector<string> binderBinSubset = MstUtils::fileToArray(binderBinSubsetPath);
                for (string binderName : binderBinSubset) {
                    Structure* seed = seedBin.getStructureNamed(binderName);
                    cout << "scoring (" << count << "/" << binderBinSubset.size() << ") " << seed->getName() << endl;

                    augmentedStructure binder(*seed);
                    scorer->setBinder(&binder);
                    scorer->defineInterfaceByPotentialContacts();

                    scorer->scoreInterface();

                    // write the score counts out to a file
                    scorer->writeContactScoresToFile();
                    scorer->writeResidueClashesToFile();
                    scorer->writeBinderScoresToFile();

                    // write out more detailed info to individual files
                    scorer->writeContactPropertyToFile(contactFilePrefix);
                    scorer->writePSSM(pssmFilePrefix);
                    seed->writePDB(pdbFilePrefix+seed->getName()+".pdb");

                    delete seed;
                    count++;
                }
            } else {
                while (seedBin.hasNext()) {
                    Structure* seed = seedBin.next();
                    cout << "scoring (" << count << "/" << totalNumSeeds << ") " << seed->getName() << endl;

                    augmentedStructure binder(*seed);
                    scorer->setBinder(&binder);
                    scorer->defineInterfaceByPotentialContacts();

                    scorer->scoreInterface();

                    // write the score counts out to a file
                    scorer->writeContactScoresToFile();
                    scorer->writeResidueClashesToFile();
                    scorer->writeBinderScoresToFile();

                    delete seed;
                    count++;
                }   
            }
        }
    }
    cout << "Done" << endl;

    return 0;
}