// ------------------------------------------------------------------
// This Code include Error State Kalman Filter and their functions
// Written By: Rasit Evduzen
// Email: rasit.evduzen@nh-tech.ai
// Version: v0.0
// Date:  27.12.2023
// Latest Update: 12.01.2023
// --------------------------------- Libraries  ---------------------------------
#include <iostream>
#include <fstream>
#include <sstream>
#include <Eigen/Dense>
#include "math_utility.hpp"
#include <vector>
// ---------------------------- Global Constants -------------------------------


// --------------------------------- Functions ---------------------------------
// -----------------------------------------------------------------------------
VectorXd deg_to_rad_v(const VectorXd &degrees) { // ! Maybe there is an Error for const referance semantics
    return degrees * PI / 180.0;
}
// -----------------------------------------------------------------------------
VectorXd rad_to_deg_v(const VectorXd &radians) {
    return radians * 180.0 / M_PI;
}
// -----------------------------------------------------------------------------
void clean_matrix_columns(MatrixXd &matrix, const vector<int>& columns, double th_min, double th_max) {
    for (int col_index = 1; col_index < columns.size(); ++col_index) {
        VectorXd column_data = matrix.col(col_index);
        double mean_value = column_data.mean();
        for (int i = 0; i < column_data.size(); ++i) {
                if (abs(column_data(i)) > th_max)
                {
                    column_data(i) = mean_value;
                    cout << "if-Data->" << column_data(i) << endl;
                    getchar();
                }
                else if (abs(column_data(i)) < (mean_value - th_min))
                {
                    column_data(i) = mean_value;
                    cout << "else if1 - Data->" << column_data(i) << endl;
                    getchar();
                } 
                else if (abs(column_data(i)) > (mean_value + th_min))
                {
                    column_data(i) = mean_value;
                    cout << "else if2 - Data->" << column_data(i) << endl;
                    getchar();
                } 
        matrix.col(col_index) = column_data;
        }
        // data_clean(column_data, th_min, th_max);
    }
}
// // -----------------------------------------------------------------------------
// void data_clean(VectorXd &data, double th_min, double th_max) {
//     double mean_value = data.mean();
//     // cout << "data_clean - selectec column " << data << endl;
//     // getchar();
//     // cout << "data_clean - column mean: " << mean_value << endl;
//     // getchar();
//     // cout << "data_clean - data size: " << data.size() << endl;
//     // getchar();

//     for (int i = 0; i < data.size(); ++i) {
//         // cout << "Data->" << data(i) << endl;
//         // getchar();
//         if (abs(data(i)) > th_max || abs(data(i)) < (mean_value - th_min) || abs(data(i)) > (mean_value + th_min))
//         {
//             data(i) = mean_value;
//             // cout << "if-Data->" << data(i) << endl;
//             // getchar();
//         }


//         // if (abs(data(i)) > th_max)
//         // {
//         //     data(i) = mean_value;
//         //     cout << "if-Data->" << data(i) << endl;
//         //     getchar();
//         // } 
//         // else if (abs(data(i)) < (mean_value - th_min))
//         // {
//         //     data(i) = mean_value;
//         //     cout << "else if1 - Data->" << data(i) << endl;
//         //     getchar();
//         // } 
//         // else if (abs(data(i)) > (mean_value + th_min))
//         // {
//         //     data(i) = mean_value;
//         //     cout << "else if1 - Data->" << data(i) << endl;
//         //     getchar();
//         // } 
//     }
// }
// -----------------------------------------------------------------------------
double calculate_mean(const VectorXd &vec) {
    // if (vec.size() == 0) {
    //     cerr << "Error Occur! Cannot calculate mean for an empty vector." << endl;
    //     return 0.0;
    // }
    return vec.mean();
}
// -----------------------------------------------------------------------------
double calculate_max(const VectorXd &vec) {
    // if (vec.size() == 0) {
    //     cerr << "Error: Cannot calculate max for an empty vector." << endl;
    //     return 0.0; 
    // }
    return vec.maxCoeff();
}
// -----------------------------------------------------------------------------
double calculate_min(const VectorXd &vec) {
    // if (vec.size() == 0) {
    //     cerr << "Error: Cannot calculate min for an empty vector." << endl;
    //     return 0.0;
    // }
    return vec.minCoeff();
}
// -----------------------------------------------------------------------------
double calculate_std_dev(const VectorXd &vec) {
    // if (vec.size() == 0) {
    //     cerr << "Error: Cannot calculate standard deviation for an empty vector." << endl;
    //     return 0.0; 
    // }
    double mean_value = vec.mean();
    double sum_squared_diffs = 0.0;
    for (int i = 0; i < vec.size(); ++i) {
        double diff = vec(i) - mean_value;
        sum_squared_diffs += diff * diff;
    }
    return sqrt(sum_squared_diffs / static_cast<double>(vec.size()));
}
// -----------------------------------------------------------------------------
double calculate_var(const VectorXd &vec) {
    // if (vec.size() == 0) {
    //     cerr << "Error: Cannot calculate variance for an empty vector." << endl;
    //     return 0.0; 
    // } 
    double mean_value = vec.mean();
    double sum_squared_diffs = 0.0;
    for (int i = 0; i < vec.size(); ++i) {
        double diff = vec(i) - mean_value;
        sum_squared_diffs += diff * diff;
    }
    return sum_squared_diffs / static_cast<double>(vec.size());
}
// -----------------------------------------------------------------------------
void accel_to_attitude(const Vector3d &accl, double &roll, double &pitch) {
    roll = atan2(-accl.y(), -accl.z());
    pitch = atan(accl.x() / sqrt(accl.y() * accl.y() + accl.z() * accl.z()));
}
// -----------------------------------------------------------------------------
void stationary_gyrocompassing(const Vector3d &W_ib_b, const Vector2d &attitude, double &yaw) {
    double roll  = attitude.x();
    double pitch = attitude.y();
    
    double w_x = W_ib_b.x();
    double w_y = W_ib_b.y();
    double w_z = W_ib_b.z();

    double sin_psi = -w_y * cos(roll)  + w_z * sin(roll);
    double cos_psi =  w_x * cos(pitch) + w_y * sin(roll) * sin(pitch) + w_z * cos(roll) * sin(pitch);

    yaw = atan2(sin_psi, cos_psi);
}
// -----------------------------------------------------------------------------
Matrix3d euler_to_ctm_XYZ(double roll, double pitch, double yaw) {
    Matrix3d ctm;
    
    double cos_roll  = cos(roll);  // rad
    double sin_roll  = sin(roll);  // rad
    double cos_pitch = cos(pitch); // rad
    double sin_pitch = sin(pitch); // rad
    double cos_yaw   = cos(yaw);   // rad
    double sin_yaw   = sin(yaw);   // rad

    ctm(0, 0) = (cos_pitch * cos_yaw);
    ctm(1, 0) = (cos_roll * sin_yaw) + (sin_roll * sin_pitch * cos_yaw);
    ctm(2, 0) = (sin_roll * sin_yaw) - (cos_roll * sin_pitch * cos_yaw);

    ctm(0, 1) = (-cos_pitch * sin_yaw);
    ctm(1, 1) = (cos_roll * cos_yaw) - (sin_roll * sin_pitch * sin_yaw);
    ctm(2, 1) = (sin_roll * cos_yaw) + (cos_roll * sin_pitch * sin_yaw);

    ctm(0, 2) = (sin_pitch); 
    ctm(1, 2) = (-sin_roll * cos_pitch);
    ctm(2, 2) = (cos_roll * cos_pitch);

    return ctm;
}
// -----------------------------------------------------------------------------
Vector3d ctm_to_euler_XYZ(const Matrix3d &ctm) {
    Vector3d euler;

    euler(1) = asin(ctm(0, 2)); // pitch 
    euler(0) = atan2(-ctm(1, 2), ctm(2, 2)); // roll
    euler(2) = atan2(-ctm(0, 1), ctm(0, 0)); // yaw 
    return euler;
}
// -----------------------------------------------------------------------------
Vector3d gravity_ned(double L_b, double h_b) {
    Vector3d gravity_vec;

    double sinsqL = sin(L_b) * sin(L_b);
    double g_0 = 9.7803253359 * (1 + 0.001931853 * sinsqL) / sqrt(1 - earth_eccen * earth_eccen * sinsqL);

    // Calculate north gravity
    gravity_vec(0) = -8.08E-9 * h_b * sin(2 * L_b);
    // East gravity is zero
    gravity_vec(1) = 0;
    // Calculate down gravity
    gravity_vec(2) = g_0 * (1 - (2 / equatorial_radius) * (1 + earth_flatten * (1 - 2 * sinsqL) + (w_ie * w_ie * equatorial_radius * equatorial_radius * polar_radius / gravity_const)) * h_b + (3 * h_b * h_b / equatorial_radius / equatorial_radius));

    return gravity_vec;
}
// -----------------------------------------------------------------------------
Matrix3d skew_symmetric(const Vector3d &vec) {
    Matrix3d skew_mat;

    skew_mat <<       0, -vec(2),  vec(1),
                 vec(2),       0, -vec(0),
                -vec(1),  vec(0),       0;

    return skew_mat;
}
// -----------------------------------------------------------------------------
Matrix3d orthonormalization(const Matrix3d &temp_matrix) {
    HouseholderQR<Matrix3d> qr(temp_matrix);
    Matrix3d matrix = qr.householderQ();
    return matrix;
}
// -----------------------------------------------------------------------------
MatrixXd state_space_matrix(const Vector3d &w_ie_w, double l_b, const Matrix3d &c_b_w, const Vector3d &f_ib_b) {
    MatrixXd ss_mat(17, 17); 
    ss_mat.block<17,17>(0,0)    = MatrixXd::Zero(17,17);
    ss_mat.block<3,3>(0,0)      = -skew_symmetric(w_ie_w);
    ss_mat.block<3,1>(0,9)      = Vector3d(0, w_ie * cos(l_b),0);
    ss_mat.block<3,1>(0,10)     = Vector3d(-w_ie * cos(l_b),0,0);
    ss_mat.block<3,3>(0,14)     = c_b_w;
    ss_mat.block<3,3>(3,0)      = -skew_symmetric(c_b_w * f_ib_b);
    ss_mat.block<3,3>(3,11)     = c_b_w;
    ss_mat.block<3,3>(6,3)      = Matrix3d::Identity();
 
    return ss_mat;
}
// -----------------------------------------------------------------------------
MatrixXd state_transition(const MatrixXd &f_mat, double Ts, int model_order) {
    int dim = f_mat.rows();
    MatrixXd Phi = MatrixXd::Identity(dim, dim);

    switch (model_order) {
        case 1:
            Phi += (f_mat * Ts);
            break;
        case 2:
            Phi += (f_mat * Ts) + (0.5 * f_mat * f_mat * Ts * Ts);
            break;
        case 3:
            Phi += (f_mat * Ts) + (0.5 * f_mat * f_mat * Ts * Ts) + ((1.0 / 6.0) * f_mat * f_mat * f_mat * Ts * Ts * Ts);
            break;
        // Add more cases for higher orders as needed
        // case 4:
        //     Phi += (f_mat * Ts) + (0.5 * f_mat * f_mat * Ts * Ts) + ((1.0 / 6.0) * f_mat * f_mat * f_mat * Ts * Ts * Ts) + ((1.0 / 24.0) * f_mat * f_mat * f_mat * f_mat * Ts * Ts * Ts * Ts);
        //     break;
        // case 5:
        //     Phi += (f_mat * Ts) + (0.5 * f_mat * f_mat * Ts * Ts) + ((1.0 / 6.0) * f_mat * f_mat * f_mat * Ts * Ts * Ts) + ((1.0 / 24.0) * f_mat * f_mat * f_mat * f_mat * Ts * Ts * Ts * Ts) + ((1.0 / 120.0) * f_mat * f_mat * f_mat * f_mat * f_mat * Ts * Ts * Ts * Ts * Ts);
        //     break;
        default:
            // cerr << "Invalid model order." << endl;
            break;
    }

    return Phi;
}
// -----------------------------------------------------------------------------
void navigation_update_wander(  double Ts, 
                                double l_b, 
                                double h_b, 
                                Matrix3d &c_b_w_minus,  
                                VectorXd &v_eb_w_minus, 
                                VectorXd &r_eb_w_minus,
                                VectorXd &f_ib_b, 
                                VectorXd &w_ib_b, 
                                double w_ang, 
                                Matrix3d &c_b_w_plus, 
                                VectorXd &v_eb_w_plus, 
                                VectorXd &r_eb_w_plus, 
                                VectorXd &w_ie_w) {

    // Earth Angular vel transform Eq [15.23]
    w_ie_w = w_ie * Vector3d(cos(w_ang) * cos(l_b), -sin(w_ang) * cos(l_b), -sin(l_b));

    // Attitude Update eq [15.18]
    Matrix3d projection_matrix;
    projection_matrix << 0,sin(l_b),-sin(w_ang) * cos(l_b),
           -sin(l_b),0,-cos(w_ang) * cos(l_b),
           sin(w_ang) * cos(l_b),cos(w_ang) * cos(l_b),0;
    c_b_w_plus = (c_b_w_minus * (Matrix3d::Identity() + skew_symmetric(w_ib_b) * Ts)) - ((w_ie * Ts) * (projection_matrix * c_b_w_minus));
    c_b_w_plus = orthonormalization(c_b_w_plus);

    // Force Transform [5.28] Update
    MatrixXd f_ib_w = 0.5 * (c_b_w_plus + c_b_w_minus) * f_ib_b;

    // Velocity Update eq[15.19]
    Vector3d grav_vec = gravity_ned(l_b, h_b); 
    v_eb_w_plus = v_eb_w_minus + (f_ib_w * Ts) + (Vector3d(0,0,grav_vec.z()) * Ts);

    // Pose Update eq[15.20]
    r_eb_w_plus = r_eb_w_minus + ((Ts / 2) * (v_eb_w_plus + v_eb_w_minus));

}
// -----------------------------------------------------------------------------
void select_columns_deg_to_rad(MatrixXd& data_matrix, const vector<int>& selected_columns, double val) {
    for (size_t i = 0; i < selected_columns.size(); ++i) {
        data_matrix.col(selected_columns[i]) = data_matrix.col(selected_columns[i]).array() * val;
    }
}
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
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
// -----------------------------------------------------------------------------
void triad_transform(const MatrixXd& source, MatrixXd& destination, const vector<int>& source_columns, const vector<int>& destination_columns) {
    // Index                   --->   [   0   1   2    3   4   5   6]
    // Raw Data Format         --->   [Time  Gx  Gy   Gz  Ax  Ay  Az]
    // Transformed Data Format --->   [Time  Gz  Gy  -GX  Az  Ay -Ax]
    
    for (size_t i = 0; i < source_columns.size(); ++i) {
        if(i == 3 || i == 6) destination.col(destination_columns[i]) = -source.col(source_columns[i]);
        else destination.col(destination_columns[i]) = source.col(source_columns[i]);
    }
}
// -----------------------------------------------------------------------------
void measurement_mean(const MatrixXd& data_matrix, VectorXd& f_ib_b_mean, VectorXd& w_ib_b_mean, double algorithm_frequency, double average_time) {
    // Data Format --->   [Time  Gz  Gy  -GX  Az  Ay -Ax]
    int row_limit = static_cast<int>(algorithm_frequency * average_time);
    
    w_ib_b_mean(0) = calculate_mean(data_matrix.topRows(row_limit).col(1));
    w_ib_b_mean(1) = calculate_mean(data_matrix.topRows(row_limit).col(2));
    w_ib_b_mean(2) = calculate_mean(data_matrix.topRows(row_limit).col(3));

    f_ib_b_mean(0) = calculate_mean(data_matrix.topRows(row_limit).col(4));
    f_ib_b_mean(1) = calculate_mean(data_matrix.topRows(row_limit).col(5));
    f_ib_b_mean(2) = calculate_mean(data_matrix.topRows(row_limit).col(6));
}
// -----------------------------------------------------------------------------
void get_offline_data(const MatrixXd& data_matrix, MatrixXd& offline_data_matrix, double algorithm_frequency, double average_time){
    int number_of_data = data_matrix.rows();
    int row_limit = (number_of_data - static_cast<int>(algorithm_frequency * average_time));
    offline_data_matrix = data_matrix.bottomRows(row_limit);
}
// -----------------------------------------------------------------------------
MatrixXd kalman_gain(const MatrixXd& P, const MatrixXd& H, const MatrixXd& R){
    return P * H.transpose() * (H * P * H.transpose() + R).inverse();  // Kalman Gain Equation!
}
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
