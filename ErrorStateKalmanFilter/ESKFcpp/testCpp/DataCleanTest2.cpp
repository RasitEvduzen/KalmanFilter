#include <iostream>
#include <Eigen/Dense>

using namespace Eigen;

void data_clean(MatrixXd &data, double th_min, double th_max) {
    for (int i = 0; i < data.rows(); ++i) {
        double mean_value = data.row(i).mean();

        for (int j = 0; j < data.cols(); ++j) {
            if (std::abs(data(i, j)) > th_max || std::abs(data(i, j)) < (mean_value - th_min) || std::abs(data(i, j)) > (mean_value + th_min)) {
                data(i, j) = mean_value;
            }
        }
    }
}

int main() {
    MatrixXd data(3, 7);
    data << 10.0, 25.0, 5.0, 30.0, 15.0, 20.0, 12.0,
            -5.0, 18.0, 7.0, 33.0, 11.0, 25.0, 8.0,
            15.0, 22.0, 6.0, 28.0, 10.0, 18.0, 14.0;

    std::cout << "Original Data:\n" << data << std::endl;

    double th_min = 5.0;
    double th_max = 25.0;
    data_clean(data, th_min, th_max);

    std::cout << "\nCleaned Data:\n" << data << std::endl;

    return 0;
}
// g++ -I C:/toolbox/eigen-3.4.0/ -o DataCleanTest2 DataCleanTest2.cpp