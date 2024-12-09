#include <iostream>
#include <Eigen/Dense>

using namespace Eigen;

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
            std::cerr << "Invalid model order." << std::endl;
            break;
    }

    return Phi;
}

int main() {

    MatrixXd F(3, 3);
    F << 1, 2, 3,
         4, 5, 6,
         7, 8, 9;

    double Ts = 0.1;
    int modelOrder = 2;

    MatrixXd result = state_transition(F, Ts, modelOrder);

    std::cout << "State Transition Matrix:\n" << result << std::endl;

    return 0;
}
// g++ -I C:/toolbox/eigen-3.4.0/ -o StateTransition StateTransition.cpp