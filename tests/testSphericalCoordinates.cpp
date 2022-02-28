#include "mstsystem.h"
#include "mstoptions.h"

int main(int argc, char *argv[]) {

    CartesianPoint testPoint(1,2,3);
    cout << "Point with coordinates: x =" << testPoint.getX() << ", y = " << testPoint.getY() << ", z = " << testPoint.getZ() << endl;

    mstreal radius = 0.0;
    mstreal polarAngle = 0.0;
    mstreal azimuthalAngle = 0.0;

    testPoint.convertToSphericalCoordinates(radius, polarAngle, azimuthalAngle);
    cout << "Point with spherical coordinates: r =" << radius << ", theta = " << polarAngle << ", psi = " << azimuthalAngle << endl;

    CartesianPoint testPoint2;
    testPoint2.setPositionBySphericalCoordinates(radius, polarAngle, azimuthalAngle);
    cout << "Cartesian to spherical to cartesian point with coordinates: x = " << testPoint2.getX() << ", y = " << testPoint2.getY() << ", z = " << testPoint2.getZ() << endl;

    return 0;
}