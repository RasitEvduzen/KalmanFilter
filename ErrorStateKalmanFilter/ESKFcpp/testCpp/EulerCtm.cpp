#include <iostream>
#include <Eigen/Dense>

using namespace Eigen;
# define PI 3.14159265358979323846   // 20 digits of pi

Matrix3d euler_to_ctm_XYZ(double roll, double pitch, double yaw) {
    Matrix3d ctm;
    
    double cos_roll  = cos(roll);  // rad
    double sin_roll  = sin(roll);  // rad
    double cos_pitch = cos(pitch); // rad
    double sin_pitch = sin(pitch); // rad
    double cos_yaw   = cos(yaw);   // rad
    double sin_yaw   = sin(yaw);   // rad

    ctm(0, 0) = (cos_pitch * cos_yaw);
    ctm(1, 0) = (cos_roll * sin_yaw) + (sin_roll * sin_pitch * cos_yaw);
    ctm(2, 0) = (sin_roll * sin_yaw) - (cos_roll * sin_pitch * cos_yaw);

    ctm(0, 1) = (-cos_pitch * sin_yaw);
    ctm(1, 1) = (cos_roll * cos_yaw) - (sin_roll * sin_pitch * sin_yaw);
    ctm(2, 1) = (sin_roll * cos_yaw) + (sin_roll * sin_pitch * sin_yaw);

    ctm(0, 2) = (sin_pitch); 
    ctm(1, 2) = (-sin_roll * cos_pitch);
    ctm(2, 2) = (cos_roll * cos_pitch);

    return ctm;
}

int main() {

    double roll  = PI/2;
    double pitch = 0;
    double yaw   = 0;

    Matrix3d ctm = euler_to_ctm_XYZ(roll, pitch, yaw);

    std::cout << "Rotation Matrix:\n" << ctm << std::endl;

    return 0;
}

// g++ -I C:/toolbox/eigen-3.4.0/ -o EulerCtm EulerCtm.cpp