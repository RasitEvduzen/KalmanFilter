#include <iostream>
#include <Eigen/Dense>

using namespace Eigen;

void data_clean(VectorXd &data, double th_min, double th_max) {
    double mean_value = data.mean();
    for (int i = 0; i < data.size(); ++i) {
        if (std::abs(data(i)) > th_max) data(i) = mean_value;
        else if (std::abs(data(i)) < (mean_value - th_min)) data(i) = mean_value;
        else if (std::abs(data(i)) > (mean_value + th_min)) data(i) = mean_value;
    }
}

int main() {
    VectorXd data(5);
    data << 8, -3, 7, 5, 2;

    double thmin = 2.0;
    double thmax = 8.0;

    std::cout << "Data Mean Value: \n" << data.mean() << "\n\n";

    std::cout << "Original data:\n" << data << "\n\n";

    data_clean(data, thmin, thmax);

    std::cout << "Cleaned data:\n" << data << "\n";

    return 0;
}
// g++ -I C:/toolbox/eigen-3.4.0/ -o DataCleanTest DataCleanTest.cpp

