#include "mstoptions.h"
#include "mstrotlib.h"
#include "mstsystem.h"

#include "generateseeds.h"
#include "residuecontact.h"
#include "residueframe.h"
#include "utilities.h"

int main(int argc, char *argv[]) {
    MstOptions op;
    op.setTitle("Loads a binary file and either A) extracts the selected structure and writes as a PDB file or B) appends structures");
    op.addOption("mode","Either 'read' or 'write'. In read-mode, will load a seedBinary file, get the selected structures, and save to PDB format. In write mode will create a new seedBinaryFile and add PDB structures (or append to an existing one) (default: read)");
    op.addOption("seedBin","Path to seed binary file. In read mode, will load this file. In write mode will either A) load the file, if one already exists or B) create a new one",true);
    op.addOption("structureName","The name of a structure. In read mode this will be extracted. In write mode this will be interpreted as a PDB file path and added to the binary file",false);
    op.addOption("structureList","The path to a file where each line is the name of a structure. In read mode each structure will be extracted from the binary file. In write mode each structure will be loaded and added to the file",false);
    op.addOption("sequenceList","The path to a file where each line is the sequence of the seed. In read mode each sequence will be threaded onto the seed before it is extracted from the binary file. In write mode each structure will be loaded, the sequence will be set, and it will be added to the file",false);
    op.addOption("targetPDB","Path to a PDB structure (read-mode only). If provided, will combine each structure from the binary file with this target structure (useful if they form a complex)",false);
    op.setOptions(argc,argv);
        
    bool read_mode = !op.isGiven("write"); //otherwise, write
    string seedBinPath = op.getString("seedBin");
    string structureName = op.getString("structureName","");
    string structureListPath = op.getString("structureList","");
    string sequenceListPath = op.getString("sequenceList","");
    string targetPDBPath = op.getString("targetPDB","");

    if ((!op.isGiven("structureName"))&(!op.isGiven("structureList"))) MstUtils::error("Must provide the name of structure, or a list of names");
    if (op.isGiven("sequenceList")&!op.isGiven("structureList")) MstUtils::error("If sequence list is provided, must also provide structure list");
    if (!read_mode && op.isGiven("targetPDB")) MstUtils::error("--targetPDB not compatible with write mode");

    if (read_mode) {
        if (!MstSys::fileExists(seedBinPath)) MstUtils::error("When in 'read' mode, must provide the path to an existing seed binary file path");

        seedBinaryFile seedBin(seedBinPath);
        seedBin.scanFilePositions();

        Structure target;
        vector<Chain*> targetChains;
        if (targetPDBPath != "") {
            target = Structure(targetPDBPath);
            for (int i = 0; i < target.chainSize(); i++) targetChains.push_back(&target.getChain(i));
        }

        if (structureName != "") {
            Structure* seed = seedBin.getStructureNamed(structureName);
            for (Chain* C : targetChains) {
                Chain* Ccopy = new Chain(*C);
                seed->appendChain(Ccopy);
            }
            seed->writePDB(seed->getName()+".pdb");
            delete seed;
        }

        if (structureListPath != "") {
            vector<string> structureList = MstUtils::fileToArray(structureListPath);
            vector<string> sequenceList = (sequenceListPath != "") ? MstUtils::fileToArray(sequenceListPath) : vector<string>();
            for (int i = 0; i < structureList.size(); i++) {
                string structureName = structureList[i];
                string sequence = (!sequenceList.empty()) ? sequenceList[i] : "";
                Structure* seed = seedBin.getStructureNamed(structureName);
                if (sequence != "") {
                    if (sequence.size() != seed->residueSize()) {
                        cout << "sequence: " << sequence.size() << " seed residues: " << seed->residueSize() << endl;
                        MstUtils::error("Number of residues in provided sequence does not match seed");
                    }
                    for (int residx = 0; residx < sequence.size(); residx++) {
                        seed->getResidue(residx).setName(SeqTools::singleToTriple(sequence.substr(residx,1)));
                    }
                }
                for (Chain* C : targetChains) {
                    Chain* Ccopy = new Chain(*C);
                    seed->appendChain(Ccopy);
                }
                seed->writePDB(seed->getName()+".pdb");
                delete seed;
            }
        }
    } else {
        seedBinaryFile* seedBin = nullptr;
        if (MstSys::fileExists(seedBinPath)) {
            cout << "Seed binary file already exists, will open file and append new structures" << endl;
            seedBin = new seedBinaryFile(seedBinPath,false,true);
        } else {
            cout << "No seed binary file is found at the provided path, will create new file" << endl;
            seedBin = new seedBinaryFile(seedBinPath,false,false);
        }

        if (structureName != "") {
            cout << "Adding one structure" << endl;
            Structure seed(structureName);
            seed.setName(MstSys::splitPath(seed.getName(),1));
            seedBin->appendStructure(&seed);
        }

        if (structureListPath != "") {
            vector<string> structureList = MstUtils::fileToArray(structureListPath);
            vector<string> sequenceList = (sequenceListPath != "") ? MstUtils::fileToArray(sequenceListPath) : vector<string>();
            if (!sequenceList.empty()) MstUtils::error("Not yet implemented");
            cout << "Adding "<< structureList.size() <<" structures" << endl;
            for (string structureName : structureList) {
                Structure seed(structureName);
                seed.setName(MstSys::splitPath(seed.getName(),1));
                seedBin->appendStructure(&seed);
            }
        }
        delete seedBin;
    }    

    cout << "Done!" << endl;
    return 0;
}