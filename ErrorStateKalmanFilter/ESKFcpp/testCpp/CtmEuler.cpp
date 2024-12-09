#include <iostream>
#include <Eigen/Dense>

using namespace Eigen;
using namespace std;

# define PI 3.14159265358979323846   // 20 digits of pi

Vector3d ctm_to_euler_XYZ(const Matrix3d &ctm) {
    Vector3d euler;

    euler(1) = asin(ctm(0, 2)); // pitch 
    euler(0) = atan2(-ctm(1, 2), ctm(2, 2)); // roll
    euler(2) = atan2(-ctm(0, 1), ctm(0, 0)); // yaw
    
    return euler;
}

int main() {

    Matrix3d ctm;
    ctm << 0.2154, -0.0215,  0.9763,
           0.9457, -0.2447, -0.2140,
           0.2435, 0.9694,  -0.0324;

    Vector3d euler = ctm_to_euler_XYZ(ctm);

    std::cout << "Euler Angles (Roll, Pitch, Yaw):\n" << euler << std::endl;

    return 0;
}
// g++ -I C:/toolbox/eigen-3.4.0/ -o CtmEuler CtmEuler.cpp