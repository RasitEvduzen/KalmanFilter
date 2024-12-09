#include <iostream>
#include <Eigen/Dense>
#include <vector>

using namespace Eigen;
using namespace std;

void select_columns_inplace(MatrixXd& data_matrix, const vector<int>& selected_columns, double scalar) {
    for (size_t i = 0; i < selected_columns.size(); ++i) {
        data_matrix.col(selected_columns[i]) = data_matrix.col(selected_columns[i]).array() * scalar;
    }
}

int main() {
    MatrixXd data_matrix(5, 7);
    data_matrix << 1, 2, 3, 4, 5, 6, 7,
                   8, 9, 10, 11, 12, 13, 14,
                   15, 16, 17, 18, 19, 20, 21,
                   22, 23, 24, 25, 26, 27, 28,
                   29, 30, 31, 32, 33, 34, 35;

    vector<int> selected_columns = {1, 3, 5};

    select_columns_inplace(data_matrix, selected_columns, 10.0);

    cout << "Processed Data:\n" << data_matrix << endl;

    return 0;
}

// g++ -I C:/toolbox/eigen-3.4.0/ -o SelectedInplace SelectedInplace.cpp 