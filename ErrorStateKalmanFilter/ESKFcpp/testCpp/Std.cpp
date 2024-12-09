#include <iostream>
#include <Eigen/Dense>
#include <cmath>

using namespace Eigen;

double calculate_std_dev(const VectorXd &vec) {
    if (vec.size() == 0) {
        std::cerr << "Error: Cannot calculate standard deviation for an empty vector." << std::endl;
        return 0.0; 
    }

    double mean_value = vec.mean();
    double sum_squared_diffs = 0.0;

    for (int i = 0; i < vec.size(); ++i) {
        double diff = vec(i) - mean_value;
        sum_squared_diffs += diff * diff;
    }

    return std::sqrt(sum_squared_diffs / static_cast<double>(vec.size()));
}

int main() {

    VectorXd inputVector(5);
    inputVector << 3.5, 1.2, 7.8, -2.3, 5.1;

    double stdDevValue = calculate_std_dev(inputVector);

    std::cout << "Standard Deviation: " << stdDevValue << std::endl;

    return 0;
}
// g++ -I C:/toolbox/eigen-3.4.0/ -o Std Std.cpp