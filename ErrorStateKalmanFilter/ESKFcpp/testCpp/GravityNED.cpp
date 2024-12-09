#include <iostream>
#include <Eigen/Dense>
#include <cmath>

using namespace Eigen;

// Constants
const double equatorial_radius = 6378137.0; // WGS84 Equatorial radius in meters
const double polar_radius = 6356752.31425; // WGS84 Polar radius in meters
const double earth_eccen = 0.0818191908425; // WGS84 eccentricity
const double earth_flatten = 1 / 298.257223563; // WGS84 earth_flattening
const double gravity_const = 3.986004418E14; // WGS84 Earth gravitational constant (m^3 s^-2)
const double w_ie = 7.292115E-5;  // Earth rotation rate (rad/s)

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

int main() {
    double l_b = 0.714242; // latitude (rad)
    double h_b = 9.9355;   // height (m)

    Vector3d gravity = gravity_ned(l_b, h_b);
    Vector3d grav_vec = gravity_ned(l_b, h_b); 
    std::cout << "Gravity (NED):\n" << gravity << std::endl;
    
    std::cout << "Gravity (NED) Vec:\n" << Vector3d(0,0,grav_vec.z()) << std::endl;

    return 0;
}
// g++ -I C:/toolbox/eigen-3.4.0/ -o GravityNED GravityNED.cpp