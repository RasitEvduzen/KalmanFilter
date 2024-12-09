#include <iostream>
#include <Eigen/Dense>

using namespace Eigen;

#define PI 3.14159265358979323846   // 20 digits of pi
#define DEG_TO_RAD (PI / 180.0)     // Scaler Degree to Radian Unit Conversion: Compile Time Constant 
#define RAD_TO_DEG (180.0 / PI)     // Scaler Radian to Degree Unit Conversion: Compile Time Constant
const double w_ie = (360. / ((((23*60)+56)*60)+4)) * DEG_TO_RAD;       // Earth Angular Velocity  [Rad/Sn]
const double equatorial_radius = 6378137;                             // Equatorial Radius Depending Standards -> WGS84    [M]   Groves Pg[56] Table 2.1
const double polar_radius = 6356752.31425;                            // Polar Radius Depending Standards -> WGS84         [M]   Groves Pg[56] Table 2.1
const double earth_flatten = 1/298.257223563;                         // Earth Flattening   [No Unit]                            Groves Pg[56] Table 2.1
const double earth_eccen = 0.0818191908425;                           // Earth Eccentricity [No Unit]                            Groves Pg[56] Table 2.1
const double gravity_const = 3.986004418E14;                          // Earth Gravity Constant   [M^3/Sn^2]

Matrix3d skew_symmetric(const Vector3d &vec) {
    Matrix3d skew_mat;

    skew_mat <<       0, -vec(2),  vec(1),
          vec(2),       0, -vec(0),
         -vec(1),  vec(0),       0;

    return skew_mat;
}
Matrix3d orthonormalization(const Matrix3d &temp_matrix) {
    HouseholderQR<Matrix3d> qr(temp_matrix);
    Matrix3d matrix = qr.householderQ();

    return matrix;
}

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

void navigation_update_wander(  double Ts, 
                                double l_b, 
                                double h_b, 
                                const Matrix3d &c_b_w_minus,  
                                const Vector3d &v_eb_w_minus, 
                                const Vector3d &r_eb_w_minus,
                                const Vector3d &f_ib_b, 
                                const Vector3d &w_ib_b, 
                                double w_ang, 
                                Matrix3d &c_b_w_plus, 
                                Vector3d &v_eb_w_plus, 
                                Vector3d &r_eb_w_plus, 
                                Vector3d &w_ie_w) {


    // Earth Angular vel transform Eq [15.23]
    w_ie_w = w_ie * Vector3d(cos(w_ang) * cos(l_b), -sin(w_ang) * cos(l_b), -sin(l_b));

    // Attitude Update eq [15.18]
    Matrix3d projection_matrix;
    projection_matrix <<       0,                sin(l_b),  -sin(w_ang) * cos(l_b),
                       -sin(l_b),                       0,  -cos(w_ang) * cos(l_b),
           sin(w_ang) * cos(l_b),   cos(w_ang) * cos(l_b),                       0;
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
int main() {

    const Matrix3d c_b_w_minus = Matrix3d::Identity();
    const Vector3d v_eb_w_minus(0.0, 0.0, 0.0);
    const Vector3d r_eb_w_minus(0.0, 0.0, 0.0);
    const Vector3d f_ib_b(1.0, 2.0, 3.0);
    const Vector3d w_ib_b(0.1, 0.2, 0.3);
    double w_ang = 0.5;

    Matrix3d c_b_w_plus;
    Vector3d v_eb_w_plus;
    Vector3d r_eb_w_plus;
    Vector3d w_ie_w_plus;

    double Ts = 0.01;
    double l_b = 0.1;
    double h_b = 100.0;

    navigation_update_wander(Ts, l_b, h_b, c_b_w_minus, v_eb_w_minus, r_eb_w_minus, f_ib_b, w_ib_b, w_ang, c_b_w_plus, v_eb_w_plus, r_eb_w_plus, w_ie_w_plus);

    std::cout << "C_b_w_plus:\n" << c_b_w_plus << std::endl;
    std::cout << "V_eb_w_plus:\n" << v_eb_w_plus << std::endl;
    std::cout << "R_eb_w_plus:\n" << r_eb_w_plus << std::endl;
    std::cout << "w_ie_w_plus: " << w_ie_w_plus << std::endl;

    //Matrix3d projection_matrix;
    //projection_matrix << 0,sin(l_b),-sin(w_ang) * cos(l_b),
    //       -sin(l_b),0,-cos(w_ang) * cos(l_b),
    //       sin(w_ang) * cos(l_b),cos(w_ang) * cos(l_b),0;

    // std::cout << "C_b_w_plus\n" << (c_bw_minus * (Matrix3d::Identity() + skew_symmetric(w_ib_b) * Ts)) - ((w_ie * Ts) * (projection_matrix * c_bw_minus)) << std::endl;

    //std::cout << projection_matrix << std::endl;


    return 0;
}
// g++ -I C:/toolbox/eigen-3.4.0/ -o NavWander  NavWander.cpp


