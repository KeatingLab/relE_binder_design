#include "mstoptions.h"
#include "mstrotlib.h"
#include "mstsystem.h"

#include "hashframes.h"
#include "generateseeds.h"
#include "residuecontact.h"
#include "residueframe.h"
#include "utilities.h"

int main(int argc, char *argv[]) {
    MstOptions op;
    op.setTitle("Loads seeds from a binary file and assigns each an address within a voxel grid");
    op.addOption("seedBin","Path to seed binary file. In read mode, will load this file. In write mode will either A) load the file, if one already exists or B) create a new one",true);
    op.addOption("voxelWidth","The side-length of the voxels in angstroms (default = 5.0 Å)",false);
    op.setOptions(argc,argv);
        
    string seedBinPath = op.getString("seedBin");
    mstreal voxelWidth = op.getReal("voxelWidth",5.0);

    if (!MstSys::fileExists(seedBinPath)) MstUtils::error("When in 'read' mode, must provide the path to an existing seed binary file path");

    seedBinaryFile seedBin(seedBinPath);
    seedBin.scanFilePositions();
    seedBin.reset();

    // Define the voxel grid
    boundingBox seedBB(0.0);
    seedBin.reset();
    while (seedBin.hasNext()) {
        Structure* seed = seedBin.next();
        CartesianPoint centroid = AtomPointerVector(seed->getAtoms()).getGeometricCenter();
        seedBB.update(centroid);
    }
    cout << "Done defining bounding box" << endl;
    seedBB.printBounds();

    // Find neighborhood address for each seed centroid and write to file
    positionHasher voxelGrid(seedBB,voxelWidth);
    fstream info_out;
    MstUtils::openFile(info_out,"seedNeighborhoodAddress.csv",fstream::out);
    info_out << "name,address" << endl;
    seedBin.reset();
    while (seedBin.hasNext()) {
        Structure* seed = seedBin.next();
        CartesianPoint centroid = AtomPointerVector(seed->getAtoms()).getGeometricCenter();
        int address = voxelGrid.hashPoint(centroid);
        info_out << seed->getName() << "," << address << endl;
    }
    info_out.close();

    cout << "Done!" << endl;
    return 0;
}