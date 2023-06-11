#include "nlohmann/json.hpp"

#include "msttypes.h"
#include "mstoptions.h"

#include "alignframes.h"
#include "generateseeds.h"
#include "scoreinterface.h"
#include "residuecontact.h"
#include "residueframe.h"

#include <chrono>

string contactPairHash(Residue *targetResA, Residue *targetResB, int distance_in_chain)
{
    return MstUtils::toString(targetResA->getResidueIndex()) + "_" + MstUtils::toString(targetResB->getResidueIndex()) + "_" + MstUtils::toString(distance_in_chain);
}

set<string> generateContactPairSet(vector<pair<Residue *, Residue *>> contacts, int binder_residue_distance = 8)
{
    // The idea: find pairs of contacts where the binder residue is close in the chain
    // For each pair of nearby contacts, then we define unique hash
    set<string> return_values;
    for (int i = 0; i < contacts.size(); i++)
    {
        pair<Residue *, Residue *> &contactA = contacts[i];
        for (int j = 0; j < contacts.size(); j++)
        {
            pair<Residue *, Residue *> &contactB = contacts[j];
            if (i == j)
                continue;
            int distance_in_binder_chain = contactA.second->getResidueIndexInChain() - contactB.second->getResidueIndexInChain();
            if (distance_in_binder_chain <= binder_residue_distance)
            {
                string hash = contactPairHash(contactA.first, contactB.first, distance_in_binder_chain);
                return_values.insert(hash);
            }
        }
    }
    return return_values;
}

int setIntersection(set<string> setA, set<string> setB)
{
    vector<string> intersection;
    std::set_intersection(setA.begin(), setA.end(), setB.begin(), setB.end(), std::back_inserter(intersection));
    return intersection.size();
}

int setUnion(set<string> setA, set<string> setB)
{
    vector<string> set_union;
    std::set_union(setA.begin(), setA.end(), setB.begin(), setB.end(), std::back_inserter(set_union));
    return set_union.size();
}

mstreal calculateJaccardSimilarity(set<string> contactPairsA, set<string> contactPairsB, int minContacts = 10)
{
    if (max(contactPairsA.size(), contactPairsB.size()) < minContacts)
        return 0.0;
    int n_intersection = setIntersection(contactPairsA, contactPairsB);
    int n_union = setUnion(contactPairsA, contactPairsB);
    mstreal sim = mstreal(n_intersection) / mstreal(n_union);
    // cout << n_intersection << " " << n_union << " " << sim << endl;
    return sim;
}

string getNameFromHeader(string header) {
    size_t pos = header.find("HEADER    ");
    if (pos == string::npos) MstUtils::error("The provided line is not a header:"+MstUtils::toString(header),"getNameFromHeader");
    return header.substr(9);
}

vector<Structure*> loadStructuresFromMultiPDB(string path, string binderChainID) {
    fstream multipdb;
    MstUtils::openFile(multipdb,path,fstream::in);
    string header_line, name;
    vector<Structure*> result;
    while (getline(multipdb,header_line)) {
        name = getNameFromHeader(header_line);
        cout << "Loading " << name << " from multiPDBfile..." << endl;
        Structure s;
        s.readPDB(multipdb);
        Chain *binderChain = s.getChainByID(binderChainID);
        if (binderChain == NULL)
            MstUtils::error("No chain matching --binderChainID found in " + MstUtils::toString(name));
        Structure* binderStructure = new Structure(*binderChain);
        binderStructure->setName(name);
        result.push_back(binderStructure);
    }
    return result;
}

vector<Structure*> loadStructuresByPath(string path, string binderChainID) {
    vector<string> binderPDBPathList = MstUtils::fileToArray(path);
    vector<Structure*> binderPDBs;
    for (string path : binderPDBPathList) {
        Structure s = Structure(path);
        Chain *binderChain = s.getChainByID(binderChainID);
        if (binderChain == NULL)
            MstUtils::error("No chain matching --binderChainID found in " + MstUtils::toString(path));
        Structure* binder = new Structure(*binderChain);
        string name = MstSys::splitPath(s.getName(), 1);
        binder->setName(name);
        binderPDBs.push_back(binder);
    }
    return binderPDBs;
}

vector<string> getBinderNames(vector<Structure*> binderStructures) {
    vector<string> result;
    for (Structure* s : binderStructures) result.push_back(s->getName());
    return result;
}

int main(int argc, char *argv[]) {
    MstOptions op;
    op.setTitle("Defines potential contacts between residues in the binder chain and residues in the target and then computes a similarity matrix.");
    op.addOption("contactData", "Path to a JSON file describing probability density of contacts", true);
    op.addOption("targetPDB", "A PDB file defining the target protein", true);
    op.addOption("structureList", "A file containing a list of PDB structures", false);
    op.addOption("multiPDB","A multiPDBfile containing multiple PDB structures separated by HEADER lines",false);
    op.addOption("binderChain", "The chain ID of the binder in the strucures (default = 0)", false);
    op.addOption("cutoff", "The cutoff used to define potential contacts (default = 0.05)", false);
    op.addOption("binderChainResidueDistance", "default = 8", false);
    op.addOption("batchID", "default = 0", false);
    op.addOption("nBatches", "default = 1", false);
    op.setOptions(argc, argv);

    string contactData = op.getString("contactData");
    string targetPDBPath = op.getString("targetPDB");
    string structureList = op.getString("structureList","");
    string multiPDB = op.getString("multiPDB","");
    string binderChainID = op.getString("binderChainID", "0");
    mstreal cutoff = op.getReal("cutoff", 0.05);
    int binderChainResidueDistance = op.getInt("binderChainResidueDistance", 8);
    int batchID = op.getInt("batchID", 0);
    int nBatches = op.getInt("nBatches", 1);

    if (nBatches < 1) MstUtils::error("nBatches < 1");
    if ((structureList=="")==(multiPDB=="")) MstUtils::error("Must provide either --structureList or --multiPDB, but not both");

    Structure targetPDB(targetPDBPath);
    vector<Chain *> binderChains;
    vector<Chain *> targetChains;

    potentialContacts potConts;
    potConts.load2DProbabilityDensities(contactData);
    potConts.setpDensityThresh(cutoff);
    potConts.setTargetResidues(targetPDB.getResidues());

    vector<Structure*> binderPDBs;
    if (structureList != "") {
        binderPDBs = loadStructuresByPath(structureList,binderChainID);
    } else {
        binderPDBs = loadStructuresFromMultiPDB(multiPDB,binderChainID);
    }
    vector<string> binderNames = getBinderNames(binderPDBs);


    vector<set<string>> contact_pair_hash_all;
    cout << "Finding contacts and generating contact pair hashes..." << endl;
    for (Structure* s : binderPDBs) {
        // cout << "Binder name: " << s.getName() << endl;
        potConts.setBinderResidues(s->getResidues());
        vector<pair<Residue *, Residue *>> contacts = potConts.getContacts();
        // cout << "Binder has " << contacts.size() << " contacts" << endl;
        // Generate contact pairs for each binder
        set<string> contact_pair_hashes = generateContactPairSet(contacts, binderChainResidueDistance);
        // cout << "Binder has " << contact_pair_hashes.size() << " contact pair hashes" << endl;
        contact_pair_hash_all.push_back(contact_pair_hashes);
    }

    // Calculate the jaccard index for each pair of binders
    if (nBatches == 1) {
        // Calculate batch boundaries
        int N_structures = binderPDBs.size();
        int N_total_calculations = N_structures * (N_structures - 1) / 2;
        cout << "In all the batches, " << N_total_calculations << " total contact pair similarity calculations" << endl;
        int batch_size = ceil(mstreal(N_total_calculations) / mstreal(nBatches));
        cout << "N batches = " << nBatches << ", so each batch size = " << batch_size << endl;
        int batch_start = (batchID)*batch_size;
        int batch_end = (batchID + 1) * batch_size;
        cout << "Batch index = " << batchID << ", so batch range = [" << batch_start << "," << batch_end << ")" << endl;

        vector<vector<mstreal>> jaccard_index_all2all;
        jaccard_index_all2all.resize(binderPDBs.size());
        int current_position = 0;
        for (int i = 0; i < contact_pair_hash_all.size(); i++)
            jaccard_index_all2all[i].resize(binderPDBs.size(), -1.0);
        for (int i = 0; i < contact_pair_hash_all.size(); i++)
        {
            if (i % 100 == 0)
            {
                cout << "i: " << i << endl;
                cout << "current position: " << current_position << endl;
            }
            for (int j = 0; j < contact_pair_hash_all.size(); j++)
            {
                if (i == j)
                    continue;
                if (i > j)
                    jaccard_index_all2all[i][j] = jaccard_index_all2all[j][i];
                // calculate new j index
                jaccard_index_all2all[i][j] = calculateJaccardSimilarity(contact_pair_hash_all[i], contact_pair_hash_all[j]);
            }
        }

        // Write data as json
        json j;
        j["binder_names"] = binderNames;
        j["jaccard_similarity"] = jaccard_index_all2all;
        fstream out;
        MstUtils::openFile(out, "binderContactJaccardSimMatrix.json", fstream::out);
        out << j.dump(4) << std::endl;
        out.close();
    }
    else {
        // Calculate batch boundaries
        int N_structures = binderPDBs.size();
        int N_total_calculations = N_structures * (N_structures - 1) / 2;
        cout << "In all the batches, " << N_total_calculations << " total contact pair similarity calculations" << endl;
        int batch_size = ceil(mstreal(N_total_calculations) / mstreal(nBatches));
        cout << "N batches = " << nBatches << ", so each batch size = " << batch_size << endl;
        int batch_start = (batchID)*batch_size;
        int batch_end = (batchID + 1) * batch_size;
        cout << "Batch index = " << batchID << ", so batch range = [" << batch_start << "," << batch_end << ")" << endl;

        map<pair<int, int>, mstreal> jaccard_index_all2all;
        int current_position = 0;
        for (int i = 0; i < contact_pair_hash_all.size(); i++) {
            for (int j = 0; j < contact_pair_hash_all.size(); j++) {
                if (current_position % 100000 == 0) {
                    cout << "i: " << i << " j: " << j << endl;
                    cout << "current position: " << current_position << endl;
                }
                if (i > j) continue;
                if (i == j) {
                    jaccard_index_all2all[pair<int, int>(i, j)] = 1.0;
                    continue;
                }
                // Note: do not increment position when on the diagonal (combinations are without replacement)
                current_position++;
                if ((current_position < batch_start) || (current_position >= batch_end)) continue;

                // calculate new j index
                jaccard_index_all2all[pair<int, int>(i, j)] = calculateJaccardSimilarity(contact_pair_hash_all[i], contact_pair_hash_all[j]);
            }
        }

        // Write data as json
        json j;
        j["binder_names"] = binderNames;
        j["jaccard_similarity"] = jaccard_index_all2all;
        fstream out;
        MstUtils::openFile(out, "binderContactJaccardSimSparse_" + MstUtils::toString(batchID) + ".json", fstream::out);
        out << j.dump(4) << std::endl;
        out.close();
    }

    for (Structure* s : binderPDBs) delete s;

    cout << "Done!" << endl;

    return 0;
}