#include <iostream>
#include <Eigen/Dense>

using namespace Eigen;

double calculate_mean(const VectorXd &vec) {
    if (vec.size() == 0) {
        std::cerr << "Error Occur! Cannot calculate mean for an empty vector." << std::endl;
        return 0.0;
    }
    return vec.mean();
}

int main() {
    VectorXd inputVector(5);
    inputVector << 3.5, 1.2, 7.8, -2.3, 5.1;

    double meanValue = calculate_mean(inputVector);

    std::cout << "Mean: " << meanValue << std::endl;

    return 0;
}
// g++ -I C:/toolbox/eigen-3.4.0/ -o Mean Mean.cpp