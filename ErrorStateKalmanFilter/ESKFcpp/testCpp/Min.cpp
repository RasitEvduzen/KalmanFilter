#include <iostream>
#include <Eigen/Dense>

using namespace Eigen;

double calculate_min(const VectorXd &vec) {
    if (vec.size() == 0) {
        std::cerr << "Error: Cannot calculate min for an empty vector." << std::endl;
        return 0.0;
    }
    return vec.minCoeff();
}



int main() {

    VectorXd inputVector(5);
    inputVector << 3.5, 1.2, 7.8, -2.3, 5.1;

    double minValue = calculate_min(inputVector);

    std::cout << "Min Value: " << minValue << std::endl;

    return 0;
}
// g++ -I C:/toolbox/eigen-3.4.0/ -o Min Min.cpp