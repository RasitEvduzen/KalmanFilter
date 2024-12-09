#include <iostream>
#include <Eigen/Dense>

using namespace Eigen;

Matrix3d skewSymmetric(const Vector3d &vec) {
    Matrix3d skew_mat;

    skew_mat <<       0, -vec(2),  vec(1),
          vec(2),       0, -vec(0),
         -vec(1),  vec(0),       0;

    return skew_mat;
}

int main() {

    Vector3d a(1.0, 2.0, 3.0);

    Matrix3d result = skewSymmetric(a);

    std::cout << "Skew-symmetric Matrix:\n" << result << std::endl;

    return 0;
}

// g++ -I C:/toolbox/eigen-3.4.0/ -o SkewMat SkewMat.cpp