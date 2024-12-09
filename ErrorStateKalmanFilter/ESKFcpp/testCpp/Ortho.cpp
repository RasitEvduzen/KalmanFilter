#include <iostream>
#include <Eigen/Dense>

using namespace Eigen;

Matrix3d orthonormalization(const Matrix3d &temp_matrix) {
    HouseholderQR<Matrix3d> qr(temp_matrix);
    Matrix3d matrix = qr.householderQ();

    return matrix;
}

int main() {
    Matrix3d tmp_matrix;
    tmp_matrix << 1.0, 2.0, 3.0,
                  4.0, 5.0, 6.0,
                  7.0, 8.0, 9.0;

    Matrix3d result = orthonormalization(tmp_matrix);
    Matrix3d t_result = result.transpose();

    std::cout << "Ortho-normalized Matrix:\n" << result << std::endl;
    std::cout << "Ortho-normalized Transpose Matrix:\n" << t_result << std::endl;
    std::cout << "Ortho-normalized Transpose Test:\n" << result*t_result << std::endl;


    return 0;
}
// g++ -I C:/toolbox/eigen-3.4.0/ -o Ortho Ortho.cpp