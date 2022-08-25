#include "msttypes.h"
#include "mstoptions.h"

#include "alignframes.h"
#include "generateseeds.h"
#include "scoreinterface.h"
#include "residuecontact.h"

#include <chrono>

int main(int argc, char *argv[]) {
    MstOptions op;
    op.setTitle("Scores a protein-protein interface by finding the prevalence of pairwise interaction geometries in the PDB");
    op.addOption("resPairDB","Path to residue pair database",true);
    op.addOption("contactData","Path to a JSON file describing probability density of contacts",true);
    op.addOption("complexPDB","A PDB file defining the complex between the protein target and designed binder. The target must have a defined sequence",false);
    op.addOption("binderChains","The binder chain ID(s) delimited by '_', e.g. '0_1'. If in --targetPDB mode, will be removed from the structure before scoring",false);
    op.addOption("targetChains","--complexPDB mode only. The protein chain(s) delimited by '_', e.g. 'A_B'",false);
    op.addOption("complexPDBList","If given, will run in complexPDB mode for multiple complexes. Argument should be followed by the path to a file where each line contains three values: [path,targetChainIDs,binderChainIDs].",false);
    op.addOption("targetPDB","A PDB file defining the structure of the protein target. Must have a defined sequence.",false);
    op.addOption("binderList","A file where each line is the path to a binder structure to be scored",false);
    op.addOption("binderBin","A seed binary file.",false);
    op.addOption("binderBinSubset","--targetPDB and --binderBin mode only. A list of structure names that will be extracted from the seed binary file and scored",false);
    op.addOption("vdwContacts","If provided, will define the interface using VDW contacts, otherwise will use CB definition. Not compatible with --targetPDB mode",false);
    op.addOption("distanceCutoff","The distance cutoff that is applied when determining whether a putative match has similar backbone atom distances (default = 0.25 Å)",false);
    // op.addOption("angleCut","The angle cutoff that is applied when a putative match has similar orientation to the query (default = 45 degrees",false);
    op.addOption("RMSDCutoff","The RMSD cutoff that is applied when confirming a putative match (default = 0.25 Å)",false);
    op.addOption("minConts","If provided, will only score binder structures with more than X contacts per residue (5 potential contacts or 2.5 vdW contacts per res)",false);
    op.addOption("verbose","");
    op.setOptions(argc,argv);

    if (op.isGiven("complexPDB")) {
        cout << "complexPDB mode" << endl;
        if (!op.isGiven("targetChains") || !op.isGiven("binderChains")) MstUtils::error("If scoring a complex, must provide '--targetChains' and '--binderChains'");
    } else if (op.isGiven("complexPDBList")) {
        cout << "complexPDBList mode" << endl;
    } else if (op.isGiven("targetPDB")) {
        cout << "targetPDB mode" << endl;
        if ((!op.isGiven("binderList")) && (!op.isGiven("binderBin"))) MstUtils::error("If providing apo target structure, must also provide '--binderList' or '--binderBin'");
        if (!op.isGiven("binderBin") && op.isGiven("binderBinSubset")) MstUtils::error("Can only extract a subset using --binderBinSubset if --binderBin is provided");
    } else {
        MstUtils::error("Must provide either --complexPDB, --complexPDBList, or --targetPDB");
    }

    string resPairDB = op.getString("resPairDB");
    string contactData = op.getString("contactData");
    string complexPDB = op.getString("complexPDB","");
    string complexPDBList = op.getString("complexPDBList","");
    string targetChainIDsString = op.getString("targetChains");
    string binderChainIDsString = op.getString("binderChains");
    string targetPDB = op.getString("targetPDB");
    string binderListPath = op.getString("binderList","");
    string seedBinPath = op.getString("binderBin","");
    string binderBinSubsetPath = op.getString("binderBinSubset","");
    bool defineVDWContacts = op.isGiven("vdwContacts");
    mstreal distanceCutoff = op.getReal("distanceCutoff",0.25);
    // mstreal angleCutoff = op.getReal("angleCut",30);
    mstreal RMSDCutoff = op.getReal("RMSDCutoff",0.25);
    mstreal minConts = op.getReal("minConts",-1);
    bool verbose = op.isGiven("verbose");

    residueBackboneBinderScorer* scorer = nullptr;

    binderScorerParams params;
    params.resPairDBPath = resPairDB;
    params.potentialContactsJSONPath = contactData;
    params.dCut = distanceCutoff;
    // params.angleCut = angleCutoff;
    params.RMSDCut = RMSDCutoff;
    params.verbose = verbose;

    if (complexPDB != "") {
        augmentedStructure complex(complexPDB,"SKIPHETERO|ALLOW ILE CD1");
        scorer = new residueBackboneBinderScorer(params,complex,binderChainIDsString,targetChainIDsString);

        if (!defineVDWContacts) {
            scorer->defineTargetBindingSiteResiduesByrSASA();
            scorer->defineInterfaceByPotentialContacts();
        }

        mstreal score = scorer->scoreBinder();
        cout << "Binder score: " << score << endl;

        scorer->writeBinderScoresToFile();
        scorer->writeContactScoresToFile();
        scorer->writeTrainingDataToFile();

    } else if (complexPDBList != "") {
        scorer = new residueBackboneBinderScorer(params);
        vector<string> lines = MstUtils::fileToArray(complexPDBList);
        for (string line : lines) {
            vector<string> val = MstUtils::split(line,",");
            if (val.size() != 3) MstUtils::error("Wrong number of values in line: "+line);
            string complexPDBPath = val[0];
            string targetChainIDsString = val[1];
            string binderChainIDsString = val[2];
            cout << complexPDBPath << " " << targetChainIDsString << " " << binderChainIDsString << endl;
            augmentedStructure complex(complexPDBPath,"SKIPHETERO|ALLOW ILE CD1");
            scorer->setComplex(complex,targetChainIDsString,binderChainIDsString);

            if (!defineVDWContacts) {
                scorer->defineTargetBindingSiteResiduesByrSASA();
                scorer->defineInterfaceByPotentialContacts();
            }

            mstreal score = scorer->scoreBinder();
            cout << "Binder score: " << score << endl;

            string listName = MstSys::splitPath(complexPDBList,1);
            scorer->writeBinderScoresToFile(true,listName);
            scorer->writeContactScoresToFile(true,listName);
            scorer->writeTrainingDataToFile(true,listName);
        }

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
        Structure target(targetS);
        scorer = new residueBackboneBinderScorer(params,target);

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
                binder.setName(MstSys::splitPath(structurePath,1));
                scorer->setBinder(&binder);
                scorer->defineInterfaceByPotentialContacts();

                mstreal score = scorer->scoreBinder();
                cout << "Binder score: " << score << endl;

                scorer->writeBinderScoresToFile();
                scorer->writeContactScoresToFile();
                scorer->writeTrainingDataToFile();

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

                    cout << "loaded " << seed->getName() << endl;

                    cout << "scoring (" << count++ << "/" << binderBinSubset.size() << ") " << seed->getName() << endl;

                    augmentedStructure binder(*seed);
                    scorer->setBinder(&binder);
                    scorer->defineInterfaceByPotentialContacts();

                    if (minConts > 0) {
                        mstreal contDensity = mstreal(scorer->countDesignableContacts())/mstreal(binder.residueSize());
                        if (contDensity < minConts) {
                            cout << "skipping binder, too few contacts" << endl;
                            delete seed;
                            continue;
                        }
                    }

                    mstreal score = scorer->scoreBinder();
                    cout << "Binder score: " << score << endl;

                    scorer->writeBinderScoresToFile();
                    scorer->writeContactScoresToFile();
                    scorer->writeTrainingDataToFile();

                    seed->writePDB(pdbFilePrefix+seed->getName()+".pdb");

                    delete seed;
                }
            } else {
                while (seedBin.hasNext()) {
                    Structure* seed = seedBin.next();
                    if (seed->residueSize() < 4) {
                        delete seed;
                        continue;
                    }
                    cout << "scoring (" << count++ << "/" << totalNumSeeds << ") " << seed->getName() << endl;

                    augmentedStructure binder(*seed);
                    scorer->setBinder(&binder);
                    scorer->defineInterfaceByPotentialContacts();

                    if (minConts > 0) {
                        mstreal contDensity = mstreal(scorer->countDesignableContacts())/mstreal(binder.residueSize());
                        if (contDensity < minConts) {
                            cout << "skipping binder, too few contacts" << endl;
                            delete seed;
                            continue;
                        }
                    }

                    mstreal score = scorer->scoreBinder();
                    cout << "Binder score: " << score << endl;

                    scorer->writeBinderScoresToFile();
                    scorer->writeContactScoresToFile();
                    scorer->writeTrainingDataToFile();

                    delete seed;
                }   
            }
        }
    }
    cout << "Done" << endl;

    return 0;
}