#include <iostream>
#include <iomanip>
#include "nlohmann/json.hpp"

#include "mstoptions.h"
#include "mstrotlib.h"
#include "mstsystem.h"

using json = nlohmann::json;

int main(int argc, char *argv[]) {
    MstOptions op;
    op.setTitle("");
    op.addOption("paramsPath", "JSON file containing the potential contact parameters",true);
    op.setOptions(argc,argv);

    string paramsPath = op.getString("paramsPath");

    // parse JSON data
    fstream fin;
    MstUtils::openFile(fin,paramsPath,fstream::in);
    json j_data = json::parse(fin);

    // see if values can be extracted
    cout << "dataPath: " << j_data["dataPath"] << endl; 
    cout << "hist0 contType: " << j_data["histList"][0]["contType"] << endl;
    cout << "hist0 AA3: " << j_data["histList"][0]["AA3"] << endl;
    cout << "hist0 var1: " << j_data["histList"][0]["var1"] << endl;
    vector<vector<int>> counts = j_data["histList"][0]["countsArray"];
    cout << "hist0 counts array with shape: " << counts.size() << "," << counts[0].size() << endl;
    cout << "histo[4][9] = " << counts[4][9] << endl;

    return 0;
}