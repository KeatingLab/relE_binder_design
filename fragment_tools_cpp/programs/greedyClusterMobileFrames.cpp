#include "msttypes.h"
#include "mstoptions.h"

#include "alignframes.h"
#include "residuecontact.h"
#include "residueframe.h"

#include <chrono>

int main(int argc, char *argv[])
{
    MstOptions op;
    op.setTitle("Greedily clusters mobile residue frames.");
    op.addOption("frameDB","Path to mobile frame database",true);
    op.addOption("posDistanceCut","Residue frames with positions within this distance can be in the same cluster",true);
    op.addOption("oriAngleCut","Residue frames with orientations within this angle can be in the same cluster",true);
    op.setOptions(argc,argv);
    return 0;
}