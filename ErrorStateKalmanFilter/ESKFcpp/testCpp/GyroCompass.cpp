#include <iostream>
#include <Eigen/Dense>

using namespace Eigen;

double stationary_gyrocompassing(const Vector3d &W_ib_b, const Vector2d &attitude) {
    double roll  = attitude.x();
    double pitch = attitude.y();
    
    double w_x = W_ib_b.x();
    double w_y = W_ib_b.y();
    double w_z = W_ib_b.z();

    double sin_psi = -w_y * cos(roll)  + w_z * sin(roll);
    double cos_psi =  w_x * cos(pitch) + w_y * sin(roll) * sin(pitch) + w_z * cos(roll) * sin(pitch);

    return atan2(sin_psi, cos_psi);
}

int main() {
    Vector3d bodyRate(0.1, -0.2, 0.3);
    Vector2d attitude(0.5, 0.2);
    std::cout << "Wx -> " << bodyRate.x() << std::endl;
    std::cout << "Wx -> " << bodyRate(0) << std::endl;

    double heading = stationary_gyrocompassing(bodyRate, attitude);

    std::cout << "Heading: " << heading << std::endl;

    return 0;
}

// g++ -I C:/toolbox/eigen-3.4.0/ -o GyroCompass GyroCompass.cpp