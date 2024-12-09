#include <iostream>
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

int main(int argc, char *argv[]) {

    Matrix <float,3,3> matrixA;
    matrixA.setZero();
    cout << matrixA << endl;

  return 0;
}
// g++ -I C:/toolbox/eigen-3.4.0/ -o main main.cpp
 
