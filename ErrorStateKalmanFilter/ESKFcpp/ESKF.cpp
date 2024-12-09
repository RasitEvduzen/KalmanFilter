// ------------------------------------------------------------------
// This Code include Error State Kalman Filter and Their Functions
// Written By: Rasit Evduzen
// Version: v0.0
// Date:  27.12.2022
// Latest Update: 03.02.2023

// --------------------------------- Libraries  ---------------------------------
#include "ESKF.hpp"
#define DEBUG_MODE 0


// ---------------------------------- Kalman Filter Initialization ----------------------------------
// Number Of State "Attitude(3x1),  Velocity(3x1), Pose(3x1), WSin(1x1), WCos(1x1), Ba(3x1), Bg(3x1)"
MatrixXd h_matrix(number_of_state, number_of_state);    // State Selection Matrix
MatrixXd f_matrix(number_of_state, number_of_state);    // State Space System Matrix
MatrixXd phi_matrix(number_of_state, number_of_state);  // State Transition Matrix
MatrixXd k_matrix(number_of_state, number_of_state);    // Kalman Gain Matrix
// ---------------------------------- P~Q~R Matrix Value ----------------------------------
// High Kalman Gain:  if R >> P => P/(P+R) ~= K ~= 0 (Algorithm belief measurement)
// Low  Kalman Gain:  if R << P => P/(P+R) ~= K ~= 1 (Algorithm belief kalman model)
MatrixXd p_matrix(number_of_state, number_of_state);    // State Covariance Matrix
double p_att = 1e-5, p_vel = 1e-4, p_pos = 1e-6, p_wander = 1e-2, p_ba = 1e-4, p_bg = 1e-12;
MatrixXd q_matrix(number_of_state, number_of_state);    // Processes Noise Covariance Matrix
double q_att = 1e-7, q_vel = 1e-7, q_pos = 1e-18, q_wander = 1e-8, q_ba = 1e-9, q_bg = 1e-19;
MatrixXd r_matrix(number_of_state, number_of_state);    // Measurement Covariance Matrix
double r_noise = 1e-9;
// ---------------------------------- Algorithm Output Parameters ----------------------------------
// ----------------------------------------------------------------
Matrix3d c_b_w_minus;       // Euler Attitude Matrix
VectorXd v_eb_w_minus(3);   // NED Frame Velocity Initial Value
VectorXd r_eb_w_minus(3);   // NED Frame Pose Initial Value
// ----------------------------------------------------------------
Matrix3d c_b_w_plus;       // Euler Attitude Matrix
VectorXd v_eb_w_plus(3);   // NED Frame Velocity Initial Value
VectorXd r_eb_w_plus(3);   // NED Frame Pose Initial Value
// ----------------------------------------------------------------
Matrix3d c_b_w_temp;       // Euler Attitude Matrix
VectorXd v_eb_w_temp(3);   // NED Frame Velocity Initial Value
VectorXd r_eb_w_temp(3);   // NED Frame Pose Initial Value
// ----------------------------------------------------------------

VectorXd w_ie_w(3);   // Wander Frame Earth Angular Veocity Transform
VectorXd b_a(3);      // Accel Dynamic Bias
VectorXd b_g(3);      // Gyro  Dynamic Bias
// ----------------------------------------------------------------

VectorXd nominal_state(number_of_state);   // Nominal State Vector (Measurement Value and Mechanization)
VectorXd error_state(number_of_state);     // Error State Vector   (ESKF)
VectorXd true_state(number_of_state);      // True State Vector = Nominal State Vector - Error State Vector
VectorXd innovation(number_of_state);      // Kalman Filter Innovation Vector

VectorXd f_ib_b(3);  // Acceleration Data
VectorXd w_ib_b(3);  // Gyro Data
VectorXd f_ib_b_mean(3); // (avarage_time * algorithm_frequency)
VectorXd w_ib_b_mean(3); // (avarage_time * algorithm_frequency)



int main(int argc, char **argv)
{
    // --------------- Set Algorithm Parameters ---------------
    h_matrix.block<17,17>(0,0)    = MatrixXd::Zero(17,17);
    h_matrix.block<3,3>(6,6)      = -Matrix3d::Identity();     // Pose State Selection
    f_matrix.block<17,17>(0,0)    = MatrixXd::Zero(17,17);
    phi_matrix.block<17,17>(0,0)  = MatrixXd::Zero(17,17);
    k_matrix.block<17,17>(0,0)    = MatrixXd::Zero(17,17);
    p_matrix.block<17,17>(0,0)    = MatrixXd::Zero(17,17);
    p_matrix.block<3,3>(0,0)      = Matrix3d::Identity() * p_att;
    p_matrix.block<3,3>(3,3)      = Matrix3d::Identity() * p_vel;
    p_matrix.block<3,3>(6,6)      = Matrix3d::Identity() * p_pos;
    p_matrix.block<2,2>(9,9)      = Matrix2d::Identity() * p_wander;
    p_matrix.block<3,3>(11,11)    = Matrix3d::Identity() * p_ba;
    p_matrix.block<3,3>(14,14)    = Matrix3d::Identity() * p_bg;
    q_matrix.block<17,17>(0,0)    = MatrixXd::Zero(17,17);
    q_matrix.block<3,3>(0,0)      = Matrix3d::Identity() * q_att;
    q_matrix.block<3,3>(3,3)      = Matrix3d::Identity() * q_vel;
    q_matrix.block<3,3>(6,6)      = Matrix3d::Identity() * q_pos;
    q_matrix.block<2,2>(9,9)      = Matrix2d::Identity() * q_wander;
    q_matrix.block<3,3>(11,11)    = Matrix3d::Identity() * q_ba;
    q_matrix.block<3,3>(14,14)    = Matrix3d::Identity() * q_bg;
    r_matrix.block<17,17>(0,0)    = MatrixXd::Zero(17,17);
    r_matrix.block<17,17>(0,0)    = MatrixXd::Identity(17,17) * r_noise;
    // ----------------------------------------------------------------
    c_b_w_minus.setZero();
    v_eb_w_minus.setZero();
    r_eb_w_minus.setZero();
    // ----------------------------------------------------------------
    c_b_w_plus.setZero();
    v_eb_w_plus.setZero();
    r_eb_w_plus.setZero();
    // ----------------------------------------------------------------
    c_b_w_temp.setZero();
    v_eb_w_temp.setZero();
    r_eb_w_temp.setZero();
    // ----------------------------------------------------------------
    w_ie_w.setZero();
    b_a.setZero();
    b_g.setZero();
    nominal_state.setZero();
    error_state.setZero();
    true_state.setZero();
    innovation.setZero();
    f_ib_b.setZero();
    w_ib_b.setZero();
    f_ib_b_mean.setZero();
    w_ib_b_mean.setZero();

    //cout << f_ib_b_mean << endl;

    // ------------------ Read Data ------------------
    // const string file_name = "data.txt"; // Data Format -> Time, Gx, Gy, Gz, Ax, Ay, Az
    const string file_name = "clean.txt"; // Data Format -> Time, Gx, Gy, Gz, Ax, Ay, Az
    int number_of_data;
    MatrixXd raw_data_matrix = read_data_file(file_name);   // Read Data From File
    vector<int> columns_to_clean = {0,1,2,3,4,5,6}; 
    cout << "Clean Data: " << raw_data_matrix << endl;
    getchar();

    clean_matrix_columns(raw_data_matrix, columns_to_clean,  th_min,  th_max); // Data Clean Error!
    cout << "Clean Data: " << raw_data_matrix << endl;
    getchar();

    // Data frame Transform  
    number_of_data = raw_data_matrix.rows();        // Number of Data
    MatrixXd data_matrix(number_of_data,7);         // Triad Transformed Data
    // Index                   --->   [   0   1   2    3   4   5   6]
    // Raw Data Format         --->   [Time  Gx  Gy   Gz  Ax  Ay  Az]
    // Transformed Data Format --->   [Time  Gz  Gy  -GX  Az  Ay -Ax]
    vector<int> source_triad_columns = {0, 3, 2, 1, 6, 5, 4}; 
    vector<int> destination_triad_columns = {0, 1, 2, 3, 4, 5, 6};
    triad_transform(raw_data_matrix, data_matrix, source_triad_columns, destination_triad_columns); // Triad Transform GO GO GO!
    vector<int> selected_columns = {1,2,3};             // Gyro Data Unit Conversion
    select_columns_deg_to_rad(data_matrix, selected_columns, DEG_TO_RAD); // Deg to rad transform
    //cout << "OK!\n" << endl;
    // cout << "---------- Raw Data ----------\n" << data_matrix << endl; // 

    // --------------------------- Initial Leveling and  Gyro Compassing -----------------------------
    measurement_mean(data_matrix, f_ib_b_mean, w_ib_b_mean, algorithm_frequency, average_time);    // Check Point! 
    // cout << "OK!\n" << endl;
    // cout << w_ib_b_mean << endl;
    // cout << "OK!\n" << endl;
    // cout << f_ib_b_mean << endl; // GO GO GO !
    accel_to_attitude(f_ib_b_mean, roll_ini, pitch_ini); // Check Point!
    // cout << "OK!\n" << endl;
    // cout << "roll-> "<< roll_ini << endl;
    // cout << "OK!\n" << endl;
    // cout << "pitch-> "<< pitch_ini << endl;
    Vector2d attitude(2);
    attitude << roll_ini, pitch_ini;
    // cout << "OK!\n" << endl;
    // cout << "Attitude-> \n"<< attitude << endl;
    stationary_gyrocompassing(w_ib_b_mean, attitude,yaw_ini); // Check Point!
    // cout << "OK!\n" << endl;
    // cout << "yaw-> \n"<< yaw_ini << endl;
    c_b_w_minus = euler_to_ctm_XYZ(roll_ini, pitch_ini, yaw_ini); // Check Point!
    Vector3d euler_initial_value;
    euler_initial_value << roll_ini, pitch_ini, yaw_ini; // Euler Angles Initial Value
    cout << "OK!\n" << endl;
    cout << "Attitude Matrix-> \n" << c_b_w_minus << endl;
    MatrixXd offline_data_matrix;
    get_offline_data(data_matrix,offline_data_matrix, algorithm_frequency, average_time); // Check Point!
    //cout << "Offline Data: \n" << offline_data_matrix << endl;
    
    MatrixXd euler_angles(offline_data_matrix.rows() + 1,3);
    // MatrixXd euler_angles(11,3);
    euler_angles.row(0) = euler_initial_value.transpose();    // Euler Angle Vector
    if (DEBUG_MODE){
        cout << "Euler Angles Initial Value:\n" << euler_angles.row(0) * RAD_TO_DEG << endl;
        cout << "Press any key to run ESKF algorithm!\n" << endl;
        getchar();
    }

    //for (size_t i = 0; i < 11; ++i)    // ESKF Calculation Loop
    for (size_t i = 0; i < (number_of_data - (static_cast<int>(algorithm_frequency * average_time))); ++i)    // ESKF Calculation Loop
    {
    if (DEBUG_MODE){ 
        cout << "loop val-> " << i << endl;
    }

    w_ib_b = offline_data_matrix.row(i).segment(1,3); // Accelerometer Data [M/Sn^2]  (Fx Fy Fz)
    f_ib_b = offline_data_matrix.row(i).rightCols(3); // Gyroscope Data [Rad/Sn]      (Wx Wy Wz)

    // ---------------------- Subtraction Estimated Bias Value ----------------------
    f_ib_b -= b_a;
    w_ib_b -= b_g;
    // ---------------------- Navigation Mechanization Wander Frame ----------------------
    if (DEBUG_MODE){
        cout << "---------- Before Mechanization ----------\n" << endl;
        cout << "w_ib_b-> \n" << w_ib_b << endl;
        cout << "f_ib_b-> \n" << f_ib_b << endl;
        cout << "c_b_w_minus-> \n" << c_b_w_minus << endl;
        cout << "v_eb_w_minus-> \n" << v_eb_w_minus << endl;
        cout << "r_eb_w_minus-> \n" << r_eb_w_minus << endl;
        cout << "Press any key to continue ESKF algorithm!\n" << endl;
        getchar();
    }

    navigation_update_wander( sampling_periode, 
                                latitude, 
                                height, 
                                c_b_w_minus,  
                                v_eb_w_minus, 
                                r_eb_w_minus,
                                f_ib_b, 
                                w_ib_b, 
                                wander_ang, 
                                c_b_w_plus, 
                                v_eb_w_plus, 
                                r_eb_w_plus, 
                                w_ie_w); // Check Point! 
    c_b_w_minus = c_b_w_plus;    // Mechanization Update
    v_eb_w_minus = v_eb_w_plus;  // Mechanization Update
    r_eb_w_minus = r_eb_w_plus;  // Mechanization Update
    
    if (DEBUG_MODE){
        cout << "---------- After Mechanization ----------\n" << endl;
        cout << "c_b_w_plus-> \n" << c_b_w_plus << endl;
        cout << "v_eb_w_plus-> \n" << v_eb_w_plus << endl;
        cout << "r_eb_w_plus-> \n" << r_eb_w_plus << endl;
        cout << "w_ie_w-> \n" << w_ie_w << endl;
        //cout << "Roll-Pitch-Yaw-> \n" << rad_to_deg_v(ctm_to_euler_XYZ(c_b_w_plus)) << endl;
        cout << "Press any key to continue ESKF algorithm!\n" << endl;
        getchar();
    }

    euler_angles.row(i + 1) = ctm_to_euler_XYZ(c_b_w_plus);
    nominal_state << ctm_to_euler_XYZ(c_b_w_plus), v_eb_w_plus, r_eb_w_plus, sin(wander_ang), cos(wander_ang), b_a, b_g;
    
    if (DEBUG_MODE){
        cout << "Nominal State Vector-> \n" << nominal_state << endl;
        cout << "Press any key to continue ESKF algorithm!\n" << endl;
        getchar();
    }

    // ------------------------- Time Update (Prediction) Stage -------------------------
    f_matrix = state_space_matrix(w_ie_w, latitude, c_b_w_plus, f_ib_b);

    if (DEBUG_MODE){
        cout << "f system matrix: \n" << f_matrix << endl;
        cout << "Press any key to continue ESKF algorithm!\n" << endl;
        getchar();
    }
    
    phi_matrix = state_transition(f_matrix, sampling_periode, model_order);

    if (DEBUG_MODE){
        cout << "State Transition Matrix:\n" << phi_matrix << endl;
        cout << "Press any key to continue ESKF algorithm!\n" << endl;
        getchar();
    }
    
    p_matrix = phi_matrix * p_matrix * phi_matrix.transpose() + q_matrix;

    if (DEBUG_MODE){
        cout << "P matrix:\n" << p_matrix << endl;
        cout << "Press any key to continue ESKF algorithm!\n" << endl;
        getchar();
    }
        
    // ------------------------- Measurement Update (Correction) Stage -------------------------
    k_matrix = kalman_gain(p_matrix, h_matrix, r_matrix); // Compute Kalman Gain! GO GO GO!
    if (DEBUG_MODE){
        cout << "K Matrix:\n" << k_matrix << endl;
        cout << "Press any key to continue ESKF algorithm!\n" << endl;
        getchar();
    }
    
    innovation.segment(6,3) = -nominal_state.segment(6,3); // Compute Innovation Vector
    if (DEBUG_MODE){
            cout << "Innovation :\n" << innovation << endl;
        cout << "Press any key to continue ESKF algorithm!\n" << endl;
        getchar();
    }
    
    error_state = error_state + k_matrix * innovation;     // Compute Error State Vector!
    if (DEBUG_MODE){
        cout << "Error State Vector:\n" << error_state << endl;
        cout << "Press any key to continue ESKF algorithm!\n" << endl;
        getchar();
    }
    
    phi_matrix = (MatrixXd::Identity(number_of_state, number_of_state) - k_matrix * h_matrix) * phi_matrix * (MatrixXd::Identity(number_of_state, number_of_state) - k_matrix * h_matrix).transpose() + k_matrix * r_matrix * k_matrix.transpose();
    
    if (DEBUG_MODE){
        cout << "P matrix:\n" << p_matrix << endl;
        cout << "Press any key to continue ESKF algorithm!\n" << endl;
        getchar();
    }
    

    // ------------------------- True State = Nominal State - Error State -------------------------
    true_state = nominal_state - error_state;
    if (DEBUG_MODE){
        cout << "True State Vector:\n" << true_state << endl;  
        cout << "Press any key to continue ESKF algorithm!\n" << endl;
        getchar();
    }
        

    // ------------------------- Close Loop Correction Phase -------------------------
    b_a = error_state.segment(11,3); // Close loop Accelerometer Dynamic Bias correction
    b_g = error_state.segment(14,3); // Close loop Gyro Dynamic Bias correction
    
    if (DEBUG_MODE){
        cout << "b_a: \n" << b_a << endl;
        cout << "b_g: \n" << b_g << endl;
        cout << "Press any key to continue ESKF algorithm!\n" << endl;
        getchar();
    }
    
    // ------------------------- True State Mechanization -------------------------
    // "Attitude(3x1),  Velocity(3x1), Pose(3x1), WSin(1x1), WCos(1x1), Ba(3x1), Bg(3x1)"
    c_b_w_temp  = euler_to_ctm_XYZ(true_state(0),true_state(1),true_state(2));    // New Mechanization Temporary Value
    c_b_w_temp = orthonormalization(c_b_w_temp);                                  // Orthonormalization
    v_eb_w_temp = true_state.segment(3,3);                                        // New Mechanization Temporary Value
    r_eb_w_temp = true_state.segment(6,3);                                        // New Mechanization Temporary Value
    wander_ang  = atan2(true_state(9),true_state(10));                            // Wander Angle [Rad]
    wander_check = ((sin(wander_ang)*sin(wander_ang)) + (cos(wander_ang)*cos(wander_ang)));   // Wander Angle Check

    navigation_update_wander( sampling_periode, 
                            latitude, 
                            height, 
                            c_b_w_temp,  
                            v_eb_w_temp, 
                            r_eb_w_temp,
                            f_ib_b, 
                            w_ib_b, 
                            wander_ang, 
                            c_b_w_plus, 
                            v_eb_w_plus, 
                            r_eb_w_plus, 
                            w_ie_w); // Check Point!

    // ------------------------- Error State Reset -------------------------
    error_state.setZero();  // Set Error State Vector to Zero

    // ------------------------- Data Logging -------------------------
    cout << "True State Euler Angles-> \n" << true_state.segment(0,3)*RAD_TO_DEG << endl;
    // cout << "Wander Angle-> \n" << wander_ang << endl;


    }
  

    // cout << "Euler Angles-> \n" << euler_angles*RAD_TO_DEG << endl;
    // cout << "True State Euler Angles-> \n" << true_state.segment(0,3)*RAD_TO_DEG << endl;


    return 0;
}


//g++ -I C:/toolbox/eigen-3.4.0/ -o ESKF ESKF.cpp math_utility.cpp
