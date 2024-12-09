#include <iostream>
#include <Eigen/Dense>

using namespace Eigen;

void accel_to_attitude(const Vector3d &accl, double &roll, double &pitch) {
    roll = atan2(-accl.y(), -accl.z());
    pitch = atan(accl.x() / sqrt(accl.y() * accl.y() + accl.z() * accl.z()));
}

int main() {
    
    Vector3d accl(0.1, -0.2, 0.3);
    std::cout << "X -> " << accl.x() << "\n";
    std::cout << "Y -> " << accl.y() << "\n";
    std::cout << "Z -> " << accl.z() << "\n";

    double roll, pitch;
    accel_to_attitude(accl, roll, pitch);

    std::cout << "Roll: " << roll << "\n";
    std::cout << "Pitch: " << pitch << std::endl;

    return 0;
}
// g++ -I C:/toolbox/eigen-3.4.0/ -o AcceltoAtt AcceltoAtt.cpp