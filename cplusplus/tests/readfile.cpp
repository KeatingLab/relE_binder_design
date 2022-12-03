#include "mstsystem.h"


int main() {
    string path = "/Users/sebastianswanson/Keating/presentations/lab_meeting/year 5/220809_labMeeting/4FXI_A_D_arg81.pdb";
    Structure S(path);
    S.writePDB("test.pdb");
    return 0;
}