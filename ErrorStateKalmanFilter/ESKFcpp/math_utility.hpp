#ifndef MATH_UTILITY_HPP // Include guard
#define MATH_UTILITY_HPP
// #pragma once // Include guard

// --------------------------------- Libraries  ---------------------------------
#include <cmath>
#include <Eigen/Dense>
#include <vector>

// --------------------------------- Name Spaces  ---------------------------------
using namespace Eigen;
using namespace std;


// --------------------------------- Macros  ---------------------------------
constexpr double PI{ 3.14159265358979323846};   // 20 digits of pi
constexpr double  DEG_TO_RAD{PI / 180.0};       // Scaler Degree to Radian Unit Conversion: Compile Time Constant 
constexpr double  RAD_TO_DEG{180.0 / PI};       // Scaler Radian to Degree Unit Conversion: Compile Time Constant

// --------------------------------- Global Constants ---------------------------------
const double w_ie = (360. / ((((23*60)+56)*60)+4)) * DEG_TO_RAD;       // Earth Angular Velocity  [Rad/Sn]
const double equatorial_radius = 6378137;                             // Equatorial Radius Depending Standards -> WGS84    [M]   Groves Pg[56] Table 2.1
const double polar_radius = 6356752.31425;                            // Polar Radius Depending Standards -> WGS84         [M]   Groves Pg[56] Table 2.1
const double earth_flatten = 1/298.257223563;                         // Earth Flattening   [No Unit]                            Groves Pg[56] Table 2.1
const double earth_eccen = 0.0818191908425;                           // Earth Eccentricity [No Unit]                            Groves Pg[56] Table 2.1
const double gravity_const = 3.986004418E14;                          // Earth Gravity Constant   [M^3/Sn^2]



// --------------------------------- Functions ---------------------------------
void clean_matrix_columns(MatrixXd &matrix, const vector<int>& columns, double th_min, double th_max);   // Data Clean via Thresholding value
VectorXd deg_to_rad_v(const VectorXd &degrees);  // Degree to Radian Unit Conversion Function For Vector Value
VectorXd rad_to_deg_v(const VectorXd &radians);  // Radian to Degree Unit Conversion Function For Vector Value
void data_clean(VectorXd &data, double th_min, double th_max);  // Data Clean Via Threasholding Method
double calculate_mean(const VectorXd &vec);      // Data mean calculate wrapper function
double calculate_max(const VectorXd &vec);       // Data max calculate wrapper function
double calculate_min(const VectorXd &vec);       // Data min calculate wrapper function
double calculate_std_dev(const VectorXd &vec);   // Data Std Dev calculate wrapper function
double calculate_var(const VectorXd &vec);       // Data Var calculate wrapper function
void accel_to_attitude(const Vector3d &accl, double &roll, double &pitch);          // Acceleration to Attitude return unit [?]
void stationary_gyrocompassing(const Vector3d &W_ib_b, const Vector2d &attitude, double &yaw); // Heading calculation      return unit [?]
Matrix3d euler_to_ctm_XYZ(double roll, double pitch, double yaw);      // Euler Angle to CTM~"DCM" matrix       input  unit [rad]
Vector3d ctm_to_euler_XYZ(const Matrix3d &ctm);                        // CTM~"DCM" to Euler Angle              output unit [rad]
Vector3d gravity_ned(double L_b, double h_b);                          // NED Frame Gravity Calculation
Matrix3d skew_symmetric(const Vector3d &vec);                           // Skew Symmetric Matrix
Matrix3d orthonormalization(const Matrix3d &temp_matrix);              // Matrix orthonormalization decomposition via QR 
MatrixXd state_transition(const MatrixXd &f_mat, double Ts, int model_order);  // State Transition Matrix "Phi" Calculation via Taylor Series 
MatrixXd state_space_matrix(const Vector3d &w_ie_w, double l_b, const Matrix3d &c_b_w, const Vector3d &f_ib_b);   // System State Space Matrix
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
                                VectorXd &w_ie_w);   // Wander Frame Mechanization Function
MatrixXd read_data_file(const string& file_name);  // Data Read From File
void select_columns_deg_to_rad(MatrixXd& data_matrix, const vector<int>& selected_columns, double scalar);  // Selected Column Operation with inplace!
void triad_transform(const MatrixXd& source, MatrixXd& destination, const vector<int>& sourceColumns, const vector<int>& destinationColumns); // Triad transform
void measurement_mean(const MatrixXd& data_matrix, VectorXd& f_ib_b_mean, VectorXd& w_ib_b_mean, double algorithm_frequency, double average_time);
void get_offline_data(const MatrixXd& data_matrix, MatrixXd& offline_data_matrix, double algorithm_frequency, double average_time);
MatrixXd kalman_gain(const MatrixXd& P, const MatrixXd& H, const MatrixXd& R);


#endif  // #ifndef MATH_UTILITY_HPP
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// ----------------------------------- UNITS -----------------------------------
// Angular Motion -> Rad, Rad/Sn, Rad/Sn^2
// Linear  Motion -> M, M/Sn, M/Sn^2
// Time           -> Sn
// Frequency      -> Hz
// 
// 
// -----------------------------------------------------------------------------
