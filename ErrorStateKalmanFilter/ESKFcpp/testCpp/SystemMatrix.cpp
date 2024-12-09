#include <iostream>
#include <Eigen/Dense>
 
using namespace Eigen;
#define PI 3.14159265358979323846   // 20 digits of pi
#define DEG_TO_RAD (PI/180)
 
//const double w_ie = (360 / ((((23*60)+56)*60)+4)) * DEG_TO_RAD;  // Earth Angular Velocity [Rad/Sn]
const double w_ie = 180 * DEG_TO_RAD;  // Earth Angular Velocity [Rad/Sn]
const int number_of_state = 17;

MatrixXd skew_symmetric(const Vector3d &a) {
    Matrix3d A;
 
    A <<     0, -a(2),  a(1),
          a(2),     0, -a(0),
         -a(1),  a(0),     0;
 
    return A;
}

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
 
int main() {
 
    Vector3d W_ie_w(0.1, 0.2, 0.3);
    //double W_ie = 0.05;
    double L_b = 0;
    Matrix3d C_b_w = Matrix3d::Identity();
    Vector3d F_ib_b (2,2,2);  // 3d vektör olcak!
 
    MatrixXd result = state_space_matrix(W_ie_w, L_b, C_b_w, F_ib_b);
 
    std::cout << "State Transition Matrix:\n" << result << std::endl;

    //std::cout << MatrixXd::Zero(number_of_state,number_of_state) << std::endl;

    // std::cout << "CBW -> \n" << C_b_w << std::endl;
    // std::cout << "F_ib_b -> \n" << F_ib_b << std::endl;
    //std::cout << "CBW*F_ib_b -> \n" << -skewSymmetric(C_b_w * F_ib_b) << std::endl;

    return 0;
}
// g++ -I C:/toolbox/eigen-3.4.0/ -o SystemMatrix SystemMatrix.cpp