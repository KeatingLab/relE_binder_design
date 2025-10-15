#include <iostream>
#include <iomanip>
#include "nlohmann/json.hpp"

#include "bridgeseeds.h"
#include "utilitiesio.h"

using json = nlohmann::json;

int main()
{
    string oldTopoDBpath = "/data1/groups/keatinglab/swans/binderDesign_relE/data/round_4_repeatNativeExtensions/4_connectFragments_relBhelix-hotspotseeds_8resConnectors_2/topoDB.json";
    topologyDB oldTopoDB;
    oldTopoDB.readDBFromJSON(oldTopoDBpath);

    string fragDBpath = "/data1/groups/keatinglab/swans/binderDesign_relE/data/round_4_repeatNativeExtensions/4_connectFragments_relBhelix-hotspotseeds_8resConnectors_2/fragmentDB.pdb";
    multiPDBFile fragDB(fragDBpath);
    fragDB.countStructures();

    string topology_name = "seed-4FXE-ARG81-relax-noHyd-48-66__0_connector-100-393-14_seed-4FXE-ARG81-relax-noHyd_E_B_E4_15654_256_267";
    // cout << oldTopoDB.data[topology_name] << endl;
    // cout << oldTopoDB.data[topology_name].size() << endl;
    // cout << oldTopoDB.data[topology_name]["connector-100-393-14"] << endl;
    // for (json::iterator it = oldTopoDB.data[topology_name].begin(); it != oldTopoDB.data[topology_name].end(); ++it) {
    //     std::cout << it.key() << " : " << it.value() << "\n";
    // }
    
    shared_ptr<fragmentTopology> test = oldTopoDB.getTopologyFromDB(topology_name,fragDB);
    test->reportTopology();
    vector<shared_ptr<topologyElement>> test2 = test->getFragments();
    cout << "vector with size " << test2.size() << endl;

    // string eh = "/data1/groups/keatinglab/swans/binderDesign_relE/analysis/round_4_2seedextensions/clusterBackbones_1seedextensions/relBhelix_1seedextensions_31resorless_cluster0.55sampled.pdb";
    // multiPDBFile seedAMulti(eh);
    // seedAMulti.loadAllStructuresSP();
    cout << "Done!" << endl;
}