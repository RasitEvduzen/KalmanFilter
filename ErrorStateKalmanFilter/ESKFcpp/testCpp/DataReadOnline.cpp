#include <iostream>
#include <fstream>
#include <sstream>
#include <Eigen/Dense>

using namespace Eigen;
using namespace std;

MatrixXd read_data_file(const string& file_name) {
    ifstream file(file_name);
    string row;
    int row_idx = 0;
    int number_of_row = 0;
    if (file.is_open()) {
        number_of_row = count(istreambuf_iterator<char>(file), istreambuf_iterator<char>(), '\n');
        file.clear();
        file.seekg(0, ios::beg);
    } else {
        cout << "File Cannot Open! " << file_name << endl;
        return MatrixXd();
    }
    MatrixXd data_matrix(number_of_row, 7);
    while (getline(file, row) && row_idx < number_of_row) {
        istringstream row_stream(row);
        string token;

        for (int col = 0; col < 7; ++col) {
            getline(row_stream, token, ',');
            data_matrix(row_idx, col) = stod(token);
        }
        row_idx++;
    }
    file.close();
    return data_matrix;
}

int main() {
    const string file_name = "data.txt";
    MatrixXd data_matrix = read_data_file(file_name);
    cout << "Data :\n" << data_matrix << endl;
    
    // if (!data_matrix.isEmpty()) {
    //     cout << "Data :\n" << data_matrix << endl;
    // }

    return 0;
}

// g++ -I C:/toolbox/eigen-3.4.0/ -o DataReadOnline DataReadOnline.cpp 