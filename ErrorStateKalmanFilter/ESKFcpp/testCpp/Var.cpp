#include <iostream>
#include <Eigen/Dense>

using namespace Eigen;

double calculate_var(const VectorXd &vec) {
    if (vec.size() == 0) {
        std::cerr << "Error: Cannot calculate variance for an empty vector." << std::endl;
        return 0.0; 
    }

    double mean_value = vec.mean();
    double sum_squared_diffs = 0.0;

    for (int i = 0; i < vec.size(); ++i) {
        double diff = vec(i) - mean_value;
        sum_squared_diffs += diff * diff;
    }

    return sum_squared_diffs / static_cast<double>(vec.size());
}

int main() {

    VectorXd inputVector(5);
    inputVector << 3.5, 1.2, 7.8, -2.3, 5.1;

    double varianceValue = calculate_var(inputVector);

    std::cout << "Variance: " << varianceValue << std::endl;

    return 0;
}
// g++ -I C:/toolbox/eigen-3.4.0/ -o Var Var.cpp