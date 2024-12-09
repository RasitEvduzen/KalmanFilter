#include <iostream>
#include <Eigen/Dense>

using namespace Eigen;

double calculate_max(const VectorXd &vec) {
    if (vec.size() == 0) {
        std::cerr << "Error: Cannot calculate max for an empty vector." << std::endl;
        return 0.0; 
    }
    return vec.maxCoeff();
}

int main() {
    VectorXd inputVector(5);
    inputVector << 3.5, 1.2, 7.8, -2.3, 5.1;

    double maxValue = calculate_max(inputVector);

    std::cout << "Max Value: " << maxValue << std::endl;

    return 0;
}
// g++ -I C:/toolbox/eigen-3.4.0/ -o Max Max.cpp